"""
This module provides functionality for accessing searched data from Proteome
Discoverer.
"""

from collections import defaultdict
from datetime import datetime
import logging
import os
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from . import (
    fetch_data, modification, paths, protein, sequence,
)


LOGGER = logging.getLogger("pyproteome.discoverer")
RE_PSP = re.compile("(\w+)\((\d+)\): ([\d\.]+)")
RE_GENE = re.compile("^>.* GN=(.+) PE=")
RE_GENE_BACKUP = re.compile("^>sp\|[\dA-Za-z]+\|([\dA-Za-z_]+) ")
RE_DESCRIPTION = re.compile(
    r"^>sp\|[\dA-Za-z]+\|[\dA-Za-z_]+ (.*?) (OS=|GN=|PE=|SV=)"
)
CONFIDENCE_MAPPING = {1: "Low", 2: "Medium", 3: "High"}


def _read_peptides(conn, pick_best_ptm=False):
    df = pd.read_sql_query(
        sql="""
        SELECT
        Peptides.PeptideID,
        Peptides.SpectrumID,
        Peptides.Sequence,
        Peptides.SearchEngineRank AS "Rank",
        Peptides.ConfidenceLevel AS "Confidence Level",
        Peptides.MissedCleavages AS "Missed Cleavages",
        PeptideScores.ScoreValue AS "Ion Score",
        SpectrumHeaders.FirstScan AS "First Scan",
        SpectrumHeaders.LastScan AS "Last Scan",
        FileInfos.FileName AS "Spectrum File",
        MassPeaks.PercentIsolationInterference AS "Isolation Interference"
        FROM Peptides
        JOIN PeptideScores
        ON Peptides.PeptideID=PeptideScores.PeptideID
        JOIN SpectrumHeaders
        ON Peptides.SpectrumID=SpectrumHeaders.SpectrumID
        JOIN FileInfos
        ON FileInfos.FileID=MassPeaks.FileID
        JOIN Masspeaks
        ON Masspeaks.MassPeakID=SpectrumHeaders.MassPeakID
        """ + ("""
        WHERE Peptides.SearchEngineRank=1
        """ if pick_best_ptm else ""),
        con=conn,
        index_col="PeptideID",
    )

    return df


def _extract_sequence(df):
    df["Sequence"] = df.apply(
        lambda row:
        sequence.extract_sequence(
            row["Proteins"],
            row["Sequence"],
        ),
        axis=1,
    )

    return df


def _extract_confidence(df):
    df["Confidence Level"] = df["Confidence Level"].apply(
        lambda x: CONFIDENCE_MAPPING[x]
    )

    return df


def _extract_spectrum_file(df):
    # "path/to/file.ext" => "file.ext"
    df["Spectrum File"] = df["Spectrum File"].apply(
        lambda x: os.path.split(x)[1]
    )

    return df


def _fix_sequence_mods(df):
    for _, row in df.iterrows():
        row["Sequence"].modifications = row["Modifications"]

    return df


def _get_proteins(df, cursor):
    prots = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        ProteinAnnotations.Description,
        Proteins.Sequence
        FROM Peptides
        JOIN PeptidesProteins
        ON Peptides.PeptideID=PeptidesProteins.PeptideID
        JOIN ProteinAnnotations
        ON ProteinAnnotations.ProteinID=PeptidesProteins.ProteinID
        JOIN Proteins
        ON Proteins.ProteinID=PeptidesProteins.ProteinID
        """,
    )

    accessions = defaultdict(list)
    genes = defaultdict(list)
    descriptions = defaultdict(list)
    sequences = defaultdict(list)

    for peptide_id, prot_string, seq in prots:
        accessions[peptide_id].append(
            fetch_data.RE_DISCOVERER_ACCESSION.match(prot_string).group(1)
        )

        gene = RE_GENE.match(prot_string)

        if not gene:
            gene = RE_GENE_BACKUP.match(prot_string)

        genes[peptide_id].append(
            gene.group(1)
        )
        descriptions[peptide_id].append(
            RE_DESCRIPTION.match(prot_string).group(1)
        )
        sequences[peptide_id].append(seq)

    # Prefetch UniProt data to speed up later queries
    # fetch_data.fetch_uniprot_data(
    #     [i for lst in accessions.values() for i in lst]
    # )

    df["Protein Descriptions"] = df.index.map(
        lambda peptide_id:
        "; ".join(descriptions[peptide_id])
    )

    df["Protein Group Accessions"] = df.index.map(
        lambda peptide_id:
        "; ".join(accessions[peptide_id])
    )

    df["Proteins"] = df.index.map(
        lambda peptide_id:
        protein.Proteins(
            proteins=tuple(
                protein.Protein(
                    accession=accession,
                    gene=gene,
                    full_sequence=seq,
                    description=desc,
                )
                for accession, gene, seq, desc in zip(
                    accessions[peptide_id],
                    genes[peptide_id],
                    sequences[peptide_id],
                    descriptions[peptide_id],
                )
            )
        )
    )

    return df


def _get_modifications(df, cursor):
    aa_mods = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        AminoAcidModifications.Abbreviation,
        PeptidesAminoAcidModifications.Position
        FROM Peptides
        JOIN PeptidesAminoAcidModifications
        ON Peptides.PeptideID=PeptidesAminoAcidModifications.PeptideID
        JOIN AminoAcidModifications
        ON PeptidesAminoAcidModifications.AminoAcidModificationID=
        AminoAcidModifications.AminoAcidModificationID
        """,
    )

    aa_mod_dict = defaultdict(list)

    for peptide_id, name, pos in aa_mods:
        try:
            sequence = df.loc[peptide_id]["Sequence"]
        except KeyError:
            continue

        mod = modification.Modification(
            rel_pos=pos,
            mod_type=name,
            nterm=False,
            cterm=False,
            sequence=sequence,
        )

        aa_mod_dict[peptide_id].append(mod)

    term_mods = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        AminoAcidModifications.Abbreviation,
        AminoAcidModifications.PositionType
        FROM Peptides
        JOIN PeptidesTerminalModifications
        ON Peptides.PeptideID=PeptidesTerminalModifications.PeptideID
        JOIN AminoAcidModifications
        ON PeptidesTerminalModifications.TerminalModificationID=
        AminoAcidModifications.AminoAcidModificationID
        """,
    )

    term_mod_dict = defaultdict(list)

    # PositionType rules taken from:
    #
    # https://github.com/compomics/thermo-msf-parser/blob/
    # 697a2fe94de2e960a9bb962d1f263dc983461999/thermo_msf_parser_API/
    # src/main/java/com/compomics/thermo_msf_parser_API/highmeminstance/
    # Parser.java#L1022
    for peptide_id, name, pos_type in term_mods:
        nterm = pos_type == 1
        pos = 0 if nterm else len(sequence)

        try:
            sequence = df.loc[peptide_id]["Sequence"]
        except KeyError:
            continue

        mod = modification.Modification(
            rel_pos=pos,
            mod_type=name,
            nterm=nterm,
            cterm=not nterm,
            sequence=sequence,
        )
        term_mod_dict[peptide_id].append(mod)

    df["Modifications"] = df.index.map(
        lambda peptide_id:
        modification.Modifications(
            mods=tuple(
                sorted(
                    term_mod_dict[peptide_id] + aa_mod_dict[peptide_id],
                    key=lambda x: (x.rel_pos, x.nterm, x.cterm),
                )
            ),
        )
    )

    return df


def _get_quantifications(df, cursor, tag_names):
    if not tag_names:
        return df

    # XXX: Bug: Peak heights do not exactly match those from Discoverer

    vals = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        ReporterIonQuanResults.QuanChannelID,
        ReporterIonQuanResults.Height
        FROM Peptides
        JOIN ReporterIonQuanResults
        ON Peptides.SpectrumID=
        ReporterIonQuanResultsSearchSpectra.SearchSpectrumID
        JOIN ReporterIonQuanResultsSearchSpectra
        ON ReporterIonQuanResultsSearchSpectra.SpectrumID=
        ReporterIonQuanResults.SpectrumID
        """,
    )

    mapping = {
        (peptide_id, channel_id): height
        for peptide_id, channel_id, height in vals
        if peptide_id in df.index
    }

    # Convert very low ion counts to nan
    for key, val in mapping.items():
        if (
            val <= 1 or
            not sequence.is_labeled(df.loc[key[0]]["Sequence"])
        ):
            mapping[key] = np.nan

    channel_ids = sorted(set(i[1] for i in mapping.keys()))

    for channel_id in channel_ids:
        tag_name = tag_names[channel_id - 1]
        df[tag_name] = df.index.map(
            lambda peptide_id:
            mapping.get((peptide_id, channel_id), np.nan)
        )

    return df


def _get_ms_data(df, cursor):
    vals = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        SpectrumHeaders.Charge,
        SpectrumHeaders.Mass,
        SpectrumHeaders.RetentionTime,
        MassPeaks.Intensity
        FROM Peptides
        JOIN SpectrumHeaders
        ON Peptides.SpectrumID=SpectrumHeaders.SpectrumID
        JOIN MassPeaks
        ON MassPeaks.MassPeakID=SpectrumHeaders.MassPeakID
        """,
    )
    mapping = {
        index: (charge, mass, rt, i)
        for index, charge, mass, rt, i in vals
    }
    charge_mapping = {
        key: {val[0]}
        for key, val in mapping.items()
    }
    mass_mapping = {
        key: {val[1]}
        for key, val in mapping.items()
    }
    rt_mapping = {
        key: {val[2]}
        for key, val in mapping.items()
    }
    i_mapping = {
        key: {val[3]}
        for key, val in mapping.items()
    }

    df["Charges"] = df.index.map(
        lambda peptide_id:
        charge_mapping.get(peptide_id, set()),
    )
    df["Masses"] = df.index.map(
        lambda peptide_id:
        mass_mapping.get(peptide_id, set()),
    )
    df["RTs"] = df.index.map(
        lambda peptide_id:
        rt_mapping.get(peptide_id, set()),
    )
    df["Intensities"] = df.index.map(
        lambda peptide_id:
        i_mapping.get(peptide_id, set()),
    )

    return df


def _get_filenames(df, cursor):
    files = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        FileInfos.PhysicalFileName
        FROM Peptides
        JOIN SpectrumHeaders
        ON Peptides.SpectrumID=SpectrumHeaders.SpectrumID
        JOIN MassPeaks
        ON MassPeaks.MassPeakID=SpectrumHeaders.MassPeakID
        JOIN FileInfos
        ON FileInfos.FileID=MassPeaks.FileID
        """,
    )
    mapping = {
        index: os.path.split(val)[-1]
        for index, val in files
    }

    df["Raw Paths"] = df.index.map(
        lambda peptide_id:
        mapping.get(peptide_id, {})
    )

    return df


def _get_q_values(df, cursor):
    fields = cursor.execute(
        """
        SELECT
        CustomDataFields.FieldID,
        CustomDataFields.DisplayName
        FROM CustomDataFields
        """,
    )
    field_ids = [
        field_id
        for field_id, name in fields
        if name in ["q-Value"]
    ]

    df["q-value"] = np.nan

    if not field_ids:
        return df

    q_vals = cursor.execute(
        """
        SELECT
        CustomDataPeptides.PeptideID,
        CustomDataPeptides.FieldValue
        FROM CustomDataPeptides
        WHERE CustomDataPeptides.FieldID IN ({})
        """.format(
            ", ".join("?" * len(field_ids))
        ),
        field_ids,
    )
    indices, vals = zip(*[
        (index, val)
        for index, val in q_vals
        if index in df.index
    ])

    df.at[indices, "q-value"] = vals

    return df


def _is_pmod(mod):
    return (
        not mod.nterm and
        not mod.cterm and
        mod.mod_type in ["Phospho"]
    )


def _reassign_mods(mods, psp_val):
    reassigned = False
    ambiguous = False

    # phophoRS example format: "T(4): 99.6; S(6): 0.4; S(10): 0.0"
    # Error messages include: "Too many isoforms"
    psp_val = [
        RE_PSP.match(i.strip())
        for i in psp_val.split(";")
    ]
    psp_val = [
        i.groups()
        for i in psp_val
        if i
    ]
    psp_val = [
        (i[0], int(i[1]), float(i[2]))
        for i in psp_val
    ]

    o_mods = [i for i in mods if not _is_pmod(i)]
    p_mods = [i for i in mods if _is_pmod(i)]
    psp_val_f = [i for i in psp_val if i[2] > 50]

    if len(p_mods) != len(psp_val_f):
        LOGGER.debug(
            "Not enough info to assign phophosite: {}".format(psp_val)
        )
        ambiguous = True
    elif set(i.rel_pos + 1 for i in p_mods) != set(i[1] for i in psp_val_f):
        p_mods = [
            modification.Modification(
                rel_pos=i[1] - 1,
                nterm=False,
                cterm=False,
                sequence=p_mods[0].sequence,
            )
            for i in psp_val_f
        ]
        reassigned = True

        mods = mods.copy()
        mods.mods = tuple(
            sorted(
                o_mods + p_mods,
                key=lambda x: (x.rel_pos, x.nterm, x.cterm)
            )
        )

    return mods, reassigned, ambiguous


def _get_phosphors(df, cursor):
    fields = cursor.execute(
        """
        SELECT
        CustomDataFields.FieldID,
        CustomDataFields.DisplayName
        FROM CustomDataFields
        """,
    )
    field_ids = [
        field_id
        for field_id, name in fields
        if name in ["phosphoRS Site Probabilities"]
    ]

    df["Ambiguous"] = False

    if not field_ids:
        return df

    psp_vals = cursor.execute(
        """
        SELECT
        CustomDataPeptides.PeptideID,
        CustomDataPeptides.FieldValue
        FROM CustomDataPeptides
        WHERE CustomDataPeptides.FieldID IN ({})
        """.format(
            ", ".join("?" * len(field_ids))
        ),
        field_ids,
    )

    changed_peptides = 0

    for pep_id, psp_val in psp_vals:
        if pep_id not in df.index:
            continue

        old_mods = df.loc[pep_id]["Modifications"]

        new_mods, reassigned, ambiguous = _reassign_mods(old_mods, psp_val)

        if reassigned:
            changed_peptides += 1

        df.at[pep_id, "Modifications"] = new_mods
        df.at[pep_id, "Ambiguous"] = ambiguous

    LOGGER.info(
        "Reassigned {} phosphosites using phosphoRS".format(changed_peptides)
    )

    return df


def read_discoverer_msf(basename, pick_best_ptm=False):
    """
    Read a Proteome Discoverer .msf file.

    Converts file contents into a pandas DataFrame similar to what one would
    get by exporting peptides to .txt directly from Discoverer.

    Parameters
    ----------
    path : str
    pick_best_ptm : bool, optional

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    msf_path = os.path.join(
        paths.MS_SEARCHED_DIR,
        basename,
    )

    LOGGER.info(
        "{}: Loading ProteomeDiscoverer peptides".format(
            os.path.splitext(basename)[0],
        )
    )
    start = datetime.now()

    with sqlite3.connect(msf_path) as conn:
        cursor = conn.cursor()

        # Get any N-terminal quantification tags
        quantification = cursor.execute(
            """
            SELECT
            ParameterValue
            FROM ProcessingNodeParameters
            WHERE ProcessingNodeParameters.ParameterName="QuantificationMethod"
            """,
        ).fetchone()

        if quantification:
            quantification = quantification[0]

            if sys.version_info.major < 3:
                quantification = quantification.encode("utf-8")

            root = ET.fromstring(quantification)
            quant_tags = root.findall(
                "MethodPart/MethodPart/Parameter[@name='TagName']",
            )
            tag_names = [i.text for i in quant_tags]
        else:
            tag_names = None

        # Read the main peptide properties
        df = _read_peptides(conn, pick_best_ptm=pick_best_ptm)

        df = _get_proteins(df, cursor)
        df = _extract_sequence(df)
        df = _extract_confidence(df)
        df = _extract_spectrum_file(df)
        df = _get_modifications(df, cursor)
        df = _get_phosphors(df, cursor)
        df = _fix_sequence_mods(df)
        df = _get_quantifications(df, cursor, tag_names)
        df = _get_q_values(df, cursor)
        df = _get_ms_data(df, cursor)
        df = _get_filenames(df, cursor)

    df["Scan Paths"] = basename

    df.reset_index(inplace=True, drop=True)

    LOGGER.info(
        "{}: Loaded {} peptides in {} hr:min:sec"
        .format(
            os.path.splitext(basename)[0],
            df.shape[0],
            str(datetime.now() - start).split('.')[0],
        )
    )

    return df
