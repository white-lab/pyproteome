"""
This module provides functionality for accessing searched data from Proteome
Discoverer.
"""

import os
import re
import sqlite3
import xml.etree.ElementTree as ET

import pandas as pd


def read_discoverer_msf(path, pep_score_cutoff=15):
    """
    Read a Proteome Discoverer .msf file.

    Converts file contents into a pandas DataFrame similar to what one would
    get by exporting peptides to .txt directly from Discoverer.

    Parameters
    ----------
    path : str
    pep_score_cutoff : int, optional

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    with sqlite3.connect(path) as conn:
        quantification = conn.cursor().execute(
            """
            SELECT
            ParameterValue
            FROM ProcessingNodeParameters
            WHERE ProcessingNodeParameters.ParameterName="QuantificationMethod"
            """,
        ).fetchone()

        if quantification:
            root = ET.fromstring(quantification[0])
            quant_tags = root.findall(
                "MethodPart/MethodPart/Parameter[@name='TagName']",
            )
            quant_titles = [i.text for i in quant_tags]
        else:
            quant_titles = None

        def _get_modifications(sequence, peptide_id):
            cursor = conn.cursor()
            aa_mods = cursor.execute(
                """
                SELECT
                AminoAcidModifications.Abbreviation,
                PeptidesAminoAcidModifications.Position
                FROM
                PeptidesAminoAcidModifications JOIN
                AminoAcidModifications
                WHERE
                PeptidesAminoAcidModifications.PeptideID=(?) AND
                PeptidesAminoAcidModifications.AminoAcidModificationID=
                AminoAcidModifications.AminoAcidModificationID
                """,
                (peptide_id,),
            )

            aa_mods = [
                "{}{}({})".format(
                    sequence[pos],
                    pos + 1,
                    name,
                )
                for name, pos in aa_mods
            ]

            term_mods = cursor.execute(
                """
                SELECT
                AminoAcidModifications.Abbreviation,
                AminoAcidModifications.PositionType
                FROM
                PeptidesTerminalModifications JOIN
                AminoAcidModifications
                WHERE
                PeptidesTerminalModifications.PeptideID=(?) AND
                PeptidesTerminalModifications.TerminalModificationID=
                AminoAcidModifications.AminoAcidModificationID
                """,
                (peptide_id,),
            )

            # PositionType rules taken from:
            #
            # https://github.com/compomics/thermo-msf-parser/blob/
            # 697a2fe94de2e960a9bb962d1f263dc983461999/thermo_msf_parser_API/
            # src/main/java/com/compomics/thermo_msf_parser_API/highmeminstance/
            # Parser.java#L1022
            term_mods = [
                "{}({})".format(
                    "N-Term" if pos_type == 1 else "C-Term",
                    name
                )
                for name, pos_type in term_mods
            ]

            mods = (
                [i for i in term_mods if i.startswith("N-Term")] +
                aa_mods +
                [i for i in term_mods if i.startswith("C-Term")]
            )

            return "; ".join(mods)

        def _get_proteins(peptide_id):
            prots = conn.cursor().execute(
                """
                SELECT
                ProteinAnnotations.Description
                FROM
                PeptidesProteins JOIN
                ProteinAnnotations
                WHERE
                PeptidesProteins.PeptideID=(?) AND
                ProteinAnnotations.ProteinID=PeptidesProteins.ProteinID
                """,
                (int(peptide_id),),
            ).fetchall()
            return "; ".join(
                re.sub(
                    r"^>sp\|[\dA-Za-z]+\|([\dA-Za-z_]+) (.*)$",
                    r"\2 - [\1]",
                    i[0],
                )
                for i in prots
            )

        def _get_quantifications(peptide_id):
            # XXX: Bug: Peak heights do not exactly match those from Discoverer
            vals = conn.cursor().execute(
                """
                SELECT
                ReporterIonQuanResults.Height
                FROM
                Peptides JOIN
                ReporterIonQuanResults JOIN
                ReporterIonQuanResultsSearchSpectra
                WHERE
                Peptides.PeptideID=(?) AND
                Peptides.SpectrumID=
                ReporterIonQuanResultsSearchSpectra.SearchSpectrumID AND
                ReporterIonQuanResultsSearchSpectra.SpectrumID=
                ReporterIonQuanResults.SpectrumID
                ORDER BY
                ReporterIonQuanResults.QuanChannelID ASC
                """,
                (peptide_id,),
            ).fetchall()

            if not vals:
                return [0] * len(quant_titles)

            # Drop any extra redundant quantifications
            vals = vals[:len(quant_titles)]
            vals = [val[0] for val in vals]

            return vals

        df = pd.read_sql_query(
            """
            SELECT
            Peptides.PeptideID,
            Peptides.Sequence,
            Peptides.ConfidenceLevel AS "Confidence Level",
            PeptideScores.ScoreValue AS "IonScore",
            SpectrumHeaders.FirstScan AS "First Scan",
            SpectrumHeaders.LastScan AS "Last Scan",
            FileInfos.FileName AS "Spectrum File"
            FROM
            Peptides JOIN
            PeptideScores JOIN
            SpectrumHeaders JOIN
            FileInfos JOIN
            Masspeaks
            WHERE
            Peptides.PeptideID=PeptideScores.PeptideID AND
            PeptideScores.ScoreValue>={} AND
            Peptides.SpectrumID=SpectrumHeaders.SpectrumID AND
            FileInfos.FileID=MassPeaks.FileID AND
            Masspeaks.MassPeakID=SpectrumHeaders.MassPeakID
            """.format(pep_score_cutoff),
            conn
        )

        # Fix >sp|P59900|EMIL3_MOUSE EMILIN-3 OS=Mus musculus GN=Emilin3 PE=2
        # SV=1
        # to EMILIN-3 OS=Mus musculus GN=Emilin3 PE=2 SV=1 - [EMIL3_MOUSE]
        df["Protein Descriptions"] = df["PeptideID"].apply(
            _get_proteins
        )

        # 1 -> "Low", 2 -> "Medium", 3 -> "High"
        confidence_mapping = {1: "Low", 2: "Medium", 3: "High"}
        df["Confidence Level"] = df["Confidence Level"].apply(
            lambda x: confidence_mapping[x]
        )

        # "path/to/file.ext" => "file.ext"
        df["Spectrum File"] = df["Spectrum File"].apply(
            lambda x: os.path.split(x)[1]
        )
        df["Modifications"] = pd.Series(
            _get_modifications(row["Sequence"], row["PeptideID"])
            for _, row in df.iterrows()
        )

        if quant_titles:
            df[quant_titles] = df.apply(
                lambda x: pd.Series(
                    _get_quantifications(x["PeptideID"])
                ),
                axis=1,
            )

    return df
