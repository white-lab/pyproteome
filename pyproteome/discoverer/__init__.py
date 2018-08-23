"""
This module provides functionality for accessing searched data from Proteome
Discoverer.
"""

from collections import defaultdict
import logging
import os
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from pyproteome import (
    pypuniprot, modification, paths, protein, sequence,
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
        SpectrumHeaders.FirstScan AS "Scan",
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
    if df.shape[0] < 1:
        return df

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
            pypuniprot.RE_DISCOVERER_ACCESSION.match(prot_string).group(1)
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

    mod_dict = defaultdict(list)

    for peptide_id, name, pos in aa_mods:
        if peptide_id not in df.index:
            continue

        mod = modification.Modification(
            rel_pos=pos,
            mod_type=name,
            nterm=False,
            cterm=False,
        )

        mod_dict[peptide_id].append(mod)

    term_mods = cursor.execute(
        """
        SELECT
        Peptides.PeptideID,
        Peptides.Sequence,
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

    # PositionType rules taken from:
    #
    # https://github.com/compomics/thermo-msf-parser/blob/
    # 697a2fe94de2e960a9bb962d1f263dc983461999/thermo_msf_parser_API/
    # src/main/java/com/compomics/thermo_msf_parser_API/highmeminstance/
    # Parser.java#L1022
    for peptide_id, pep_seq, name, pos_type in term_mods:
        if peptide_id not in df.index:
            continue

        nterm = pos_type == 1
        pos = 0 if nterm else len(pep_seq)

        mod = modification.Modification(
            rel_pos=pos,
            mod_type=name,
            nterm=nterm,
            cterm=not nterm,
        )
        mod_dict[peptide_id].append(mod)

    mod_dict = {
        key: _sort_mods(val)
        for key, val in mod_dict.items()
    }

    def _get_mods(row):
        peptide_id = row.name

        mods = modification.Modifications(
            mods=mod_dict.get(peptide_id, tuple()),
        )

        for mod in mods.mods:
            assert mod.sequence is None
            mod.sequence = row["Sequence"]

        row["Sequence"].modifications = mods

        return mods

    df["Modifications"] = df.apply(_get_mods, axis=1)

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

    channel_ids = sorted(set(i[1] for i in mapping.keys()))

    col_names = [tag_names[channel_id - 1] for channel_id in channel_ids]

    def _get_quants(row):
        peptide_id = row.name

        vals = [
            mapping.get((peptide_id, channel_id), np.nan)
            for channel_id in channel_ids
        ]

        # Convert very low ion counts and unlabeled peptides to nan
        vals = [
            np.nan if val <= 1 else val
            for val in vals
        ]

        if not np.isnan(vals).all() and not row["Sequence"].is_labeled:
            vals = [np.nan] * len(vals)

        return pd.Series(vals, index=col_names)

    df[col_names] = df.apply(
        _get_quants,
        axis=1,
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


def _set_defaults(df):
    df["Validated"] = False
    df["Fold Change"] = np.nan
    df["p-value"] = np.nan

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


def _sort_mods(mods):
    return tuple(
        sorted(
            mods,
            key=lambda x: (x.rel_pos, x.nterm, x.cterm, x.mod_type),
        )
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
                mod_type="Phospho",
                nterm=False,
                cterm=False,
                sequence=p_mods[0].sequence,
            )
            for i in psp_val_f
        ]
        reassigned = True

        mods = modification.Modifications(
            mods=_sort_mods(o_mods + p_mods),
        )

        for mod in mods.mods:
            mod.sequence.modifications = mods

    return mods, reassigned, ambiguous


def _get_phosphors(df, cursor, name=None):
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
        "{}: -- Reassigned {} phosphosites using phosphoRS".format(
            name,
            changed_peptides,
        )
    )

    return df


def _get_species(cursor):
    fields = cursor.execute(
        """
        SELECT
        ProcessingNodeParameters.ParameterValue
        FROM ProcessingNodeParameters
        WHERE ProcessingNodeParameters.ParameterName="Taxonomy"
        """,
    )
    species = set()

    for name, in fields:
        # ". . . . . . . . . . . . . . . . Homo sapiens (human)"
        # => "Homo sapiens"
        species.add(" ".join(name.strip(". ").split(" ")[:2]))

    return species


def _update_label_names(cursor):
    fields = cursor.execute(
        """
        SELECT
        AminoAcidModifications.Abbreviation,
        AminoAcidModifications.ModificationName,
        AminoAcids.AminoAcidName,
        AminoAcids.OneLetterCode
        FROM AminoAcidModifications
        JOIN AminoAcidModificationsAminoAcids
        ON AminoAcidModifications.AminoAcidModificationID=
        AminoAcidModificationsAminoAcids.AminoAcidModificationID
        JOIN AminoAcids
        ON AminoAcids.AminoAcidID=
        AminoAcidModificationsAminoAcids.AminoAcidID
        WHERE AminoAcidModifications.isActive
        """,
    )

    for abbrev, mod_name, aa_name, letter in fields:
        if any(
            i in abbrev or i in mod_name
            for i in modification.LABEL_NAME_TARGETS
        ):
            if aa_name in ["N-Terminus"]:
                letter = "N-term"
            elif aa_name in ["C-Terminus"]:
                letter = "C-term"

            modification.LABEL_NAMES[abbrev].add(letter)


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
    df : :class:`pandas.DataFrame`
    """
    msf_path = os.path.join(
        paths.MS_SEARCHED_DIR,
        basename,
    )

    if not os.path.exists(msf_path):
        raise Exception("Search database does not exist: {}".format(msf_path))

    name = os.path.splitext(basename)[0]

    LOGGER.info(
        "{}: Loading ProteomeDiscoverer peptides...".format(name)
    )

    with sqlite3.connect(msf_path) as conn:
        cursor = conn.cursor()

        _update_label_names(cursor)

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
        df = _get_phosphors(df, cursor, name=name)
        df = _get_q_values(df, cursor)
        df = _get_ms_data(df, cursor)
        df = _get_filenames(df, cursor)
        df = _get_quantifications(df, cursor, tag_names)

        df = _set_defaults(df)

        species = _get_species(cursor)

    df["Scan Paths"] = basename

    df.reset_index(inplace=True, drop=True)

    LOGGER.info(
        "{}: Loaded {} peptides"
        .format(
            os.path.splitext(basename)[0],
            df.shape[0],
        )
    )

    return df, species
