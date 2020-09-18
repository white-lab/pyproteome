'''
This module provides functionality for accessing searched data from Proteome
Discoverer.
'''

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
    data_sets, paths, pypuniprot,
)


LOGGER = logging.getLogger('pyproteome.discoverer')
RE_PSP = re.compile(r'(\w+)\((\d+)\): ([\d\.]+)')
RE_GENE = re.compile(r'^>.* GN=([A-Za-z0-9_\-]+)( \w+)?')
RE_GENE_BACKUP = re.compile(r'^>(gi\|[\dA-Za-z]+\|)?sp\|[\dA-Za-z\-]+\|([\dA-Za-z_]+) ')
RE_DESCRIPTION = re.compile(
    r'^>(gi\|[\dA-Za-z]+\|)?sp\|[\dA-Za-z\-]+\|[\dA-Za-z_]+ (.*?) (OS=|GN=|PE=|SV=)'
)
CONFIDENCE_MAPPING = {1: 'Low', 2: 'Medium', 3: 'High'}


def _get_pd_version(cursor):
    query = cursor.execute(
        '''
        SELECT
        SoftwareVersion

        FROM SchemaInfo

        WHERE Kind=='Result'
        ''',
    )
    return tuple([int(i) for i in next(query)[0].split('.')])


def _update_label_names(cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        fields = cursor.execute(
            '''
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
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        fields = cursor.execute(
            '''
            SELECT
            FoundModifications.Abbreviation,
            FoundModifications.Name,
            AminoAcids.Name,
            AminoAcids.OneLetterCode

            FROM FoundModifications

            JOIN FoundModificationsAminoAcids
            ON FoundModifications.ModificationID=
            FoundModificationsAminoAcids.FoundModificationsModificationID

            JOIN AminoAcids
            ON AminoAcids.AminoAcidID=
            FoundModificationsAminoAcids.AminoAcidsAminoAcidID
            ''',
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    for abbrev, mod_name, aa_name, letter in fields:
        if any(
            i in abbrev or i in mod_name
            for i in data_sets.modification.LABEL_NAME_TARGETS
        ):
            if aa_name in ['N-Terminus']:
                letter = 'N-term'
            elif aa_name in ['C-Terminus']:
                letter = 'C-term'

            data_sets.modification.LABEL_NAMES[abbrev].add(letter)


def _get_quant_tags(cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        quantification = cursor.execute(
            '''
            SELECT
            ParameterValue

            FROM ProcessingNodeParameters

            WHERE ProcessingNodeParameters.ParameterName='QuantificationMethod'
            ''',
        ).fetchone()

        if quantification:
            quantification = quantification[0]

            if sys.version_info.major < 3:
                quantification = quantification.encode('utf-8')

            root = ET.fromstring(quantification)
            quant_tags = root.findall(
                'MethodPart/MethodPart/Parameter[@name=\'TagName\']',
            )
            tag_names = [i.text for i in quant_tags]
        else:
            tag_names = None
    elif pd_version[:2] in [(2, 2)]:
        quantification = cursor.execute(
            '''
            SELECT
            AnalysisDefinitionXML

            FROM AnalysisDefinition
            ''',
        ).fetchone()

        if quantification:
            quantification = quantification[0]

            if sys.version_info.major < 3:
                quantification = quantification.encode('utf-8')

            root = ET.fromstring(quantification)
            quant_tags = root.findall(
                'StudyDefinition/QuanMethods/'
                'QuanMethod/QuanChannels/QuanChannel',
            )
            tag_names = [i.get('Name') for i in quant_tags]
        else:
            tag_names = None
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    return tag_names


def _read_peptides(conn, pd_version, pick_best_psm=False):
    if pd_version[:2] in [(1, 4)]:
        sql_query = '''
        SELECT
        Peptides.PeptideID,
        Peptides.SpectrumID,
        Peptides.Sequence,
        Peptides.SearchEngineRank AS 'Rank',
        Peptides.ConfidenceLevel AS 'Confidence Level',
        Peptides.MissedCleavages AS 'Missed Cleavages',
        PeptideScores.ScoreValue AS 'Ion Score',
        SpectrumHeaders.FirstScan AS 'Scan',
        SpectrumHeaders.LastScan AS 'Last Scan',
        FileInfos.FileName AS 'Spectrum File',
        MassPeaks.PercentIsolationInterference AS 'Isolation Interference'

        FROM Peptides

        JOIN PeptideScores
        ON Peptides.PeptideID=PeptideScores.PeptideID

        JOIN SpectrumHeaders
        ON Peptides.SpectrumID=SpectrumHeaders.SpectrumID

        JOIN FileInfos
        ON FileInfos.FileID=MassPeaks.FileID

        JOIN Masspeaks
        ON Masspeaks.MassPeakID=SpectrumHeaders.MassPeakID
        ''' + (
            '''
            WHERE Peptides.SearchEngineRank=1
            '''
            if pick_best_psm else
            ''
        )
    elif pd_version[:2] in [(2, 2)]:
        sql_query = '''
        SELECT
        TargetPsms.PeptideID,
        TargetPsmsMSnSpectrumInfo.MSnSpectrumInfoSpectrumID,
        TargetPsms.Sequence,
        TargetPsms.SearchEngineRank AS 'Rank',
        TargetPsms.MissedCleavages AS 'Missed Cleavages',
        TargetPsms.IonsScore AS 'Ion Score',
        TargetPsms.FirstScan AS 'Scan',
        TargetPsms.LastScan AS 'Last Scan',
        WorkflowInputFiles.FileName AS 'Spectrum File',
        TargetPsms.PercentIsolationInterference AS 'Isolation Interference'

        FROM TargetPsms

        JOIN TargetPsmsMSnSpectrumInfo
        ON TargetPsmsMSnSpectrumInfo.TargetPsmsPeptideID=TargetPsms.PeptideID

        JOIN WorkflowInputFiles
        ON WorkflowInputFiles.FileID=TargetPsms.SpectrumFileID
        ''' + (
            '''
            WHERE TargetPsms.SearchEngineRank=1
            '''
            if pick_best_psm else
            ''
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    df = pd.read_sql_query(
        sql=sql_query,
        con=conn,
        index_col='PeptideID',
    )

    if 'Confidence Level' not in df.columns:
        # XXX: Hacked on!
        df['Confidence Level'] = 3

    return df


def _extract_sequence(df):
    if df.shape[0] < 1:
        return df

    df['Sequence'] = df.apply(
        lambda row:
        data_sets.extract_sequence(
            row['Proteins'],
            row['Sequence'],
        ),
        axis=1,
    )

    return df


def _extract_confidence(df):
    df['Confidence Level'] = df['Confidence Level'].apply(
        lambda x: CONFIDENCE_MAPPING[x]
    )

    return df


def _extract_spectrum_file(df):
    # 'path/to/file.ext' => 'file.ext'
    df['Spectrum File'] = df['Spectrum File'].apply(
        lambda x: os.path.split(x)[1]
    )

    return df


def _get_proteins(df, cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        prots = cursor.execute(
            '''
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
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        prots = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            TargetProteins.FastaTitleLines,
            TargetProteins.Sequence

            FROM TargetPsms

            JOIN TargetProteinsTargetPsms
            ON TargetPsms.PeptideID=
            TargetProteinsTargetPsms.TargetPsmsPeptideID

            JOIN TargetProteins
            ON TargetProteins.UniqueSequenceID=
            TargetProteinsTargetPsms.TargetProteinsUniqueSequenceID
            ''',
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    accessions = defaultdict(list)
    genes = defaultdict(list)
    descriptions = defaultdict(list)
    sequences = defaultdict(list)

    for peptide_id, prot_string, seq in prots:
        for fasta_line in prot_string.split('\n'):
            try:
                accessions[peptide_id].append(
                    pypuniprot.RE_DISCOVERER_ACCESSION.match(fasta_line).group(2)
                )
            except:
                print(fasta_line)
                raise

            gene = RE_GENE.match(prot_string)

            if gene:
                gene = gene.group(1)
            else:
                gene = RE_GENE_BACKUP.match(prot_string).group(2)

            genes[peptide_id].append(
                gene
            )
            descriptions[peptide_id].append(
                RE_DESCRIPTION.match(prot_string).group(2)
            )
            sequences[peptide_id].append(seq)

    df['Protein Descriptions'] = df.index.map(
        lambda peptide_id:
        '; '.join(descriptions[peptide_id])
    )

    df['Protein Group Accessions'] = df.index.map(
        lambda peptide_id:
        '; '.join(accessions[peptide_id])
    )

    df['Proteins'] = df.index.map(
        lambda peptide_id:
        data_sets.Proteins(
            proteins=tuple(
                data_sets.Protein(
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


def _get_modifications(df, cursor, pd_version):
    mod_dict = defaultdict(list)

    if pd_version[:2] in [(1, 4)]:
        aa_mods = cursor.execute(
            '''
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
            ''',
        )

        for peptide_id, name, pos in aa_mods:
            if peptide_id not in df.index:
                continue

            mod = data_sets.Modification(
                rel_pos=pos,
                mod_type=name,
                nterm=False,
                cterm=False,
            )

            mod_dict[peptide_id].append(mod)

        term_mods = cursor.execute(
            '''
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
            ''',
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

            mod = data_sets.Modification(
                rel_pos=pos,
                mod_type=name,
                nterm=nterm,
                cterm=not nterm,
            )
            mod_dict[peptide_id].append(mod)
    elif pd_version[:2] in [(2, 2)]:
        aa_mods = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            FoundModifications.Abbreviation,
            TargetPsmsFoundModifications.Position

            FROM TargetPsms

            JOIN TargetPsmsFoundModifications
            ON
            TargetPsmsFoundModifications.TargetPsmsPeptideID=TargetPsms.PeptideID

            JOIN FoundModifications
            ON
            TargetPsmsFoundModifications.FoundModificationsModificationID=
            FoundModifications.ModificationID

            WHERE
            FoundModifications.PositionType NOT IN (1, 2)
            ''',
        )

        for peptide_id, name, pos in aa_mods:
            if peptide_id not in df.index:
                continue

            pos -= 1

            mod = data_sets.Modification(
                rel_pos=pos,
                mod_type=name,
                nterm=False,
                cterm=False,
            )

            mod_dict[peptide_id].append(mod)

        term_mods = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            TargetPsms.Sequence,
            FoundModifications.Abbreviation,
            FoundModifications.PositionType

            FROM TargetPsms

            JOIN TargetPsmsFoundModifications
            ON
            TargetPsmsFoundModifications.TargetPsmsPeptideID=TargetPsms.PeptideID

            JOIN FoundModifications
            ON
            TargetPsmsFoundModifications.FoundModificationsModificationID=
            FoundModifications.ModificationID

            WHERE
            FoundModifications.PositionType IN (1, 2)
            ''',
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

            mod = data_sets.Modification(
                rel_pos=pos,
                mod_type=name,
                nterm=nterm,
                cterm=not nterm,
            )
            mod_dict[peptide_id].append(mod)
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    mod_dict = {
        key: _sort_mods(val)
        for key, val in mod_dict.items()
    }

    def _get_mods(row):
        peptide_id = row.name

        mods = data_sets.Modifications(
            mods=mod_dict.get(peptide_id, tuple()),
        )

        for mod in mods.mods:
            assert mod.sequence is None
            mod.sequence = row['Sequence']

        row['Sequence'].modifications = mods

        return mods

    df['Modifications'] = df.apply(_get_mods, axis=1)

    return df


def _get_quantifications(df, cursor, pd_version, tag_names):
    if not tag_names:
        return df

    # XXX: Bug: Peak heights do not exactly match those from Discoverer

    if pd_version[:2] in [(1, 4)]:
        vals = cursor.execute(
            '''
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
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        vals = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            ReporterQuanResults.QuanChannelID,
            ReporterQuanResults.Height

            FROM ReporterQuanResults

            INNER JOIN QuanSpectrumInfoTargetPsms
            ON QuanSpectrumInfoTargetPsms.QuanSpectrumInfoSpectrumID=
            ReporterQuanResults.QuanResultID

            INNER JOIN QuanSpectrumInfo
            ON QuanSpectrumInfo.SpectrumID=ReporterQuanResults.QuanResultID

            INNER JOIN TargetPsms
            ON TargetPsms.PeptideID=
            QuanSpectrumInfoTargetPsms.TargetPsmsPeptideID
            ''',
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
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

        if not row['Sequence'].is_labeled:
            vals = [np.nan] * len(vals)

        return pd.Series(vals, index=col_names)

    df[col_names] = df.apply(
        _get_quants,
        axis=1,
    )

    return df


def _get_ms_data(df, cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        vals = cursor.execute(
            '''
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
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        vals = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            TargetPsms.Charge,
            TargetPsms.Mass,
            TargetPsms.RetentionTime,
            TargetPsms.Intensity

            FROM TargetPsms
            ''',
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
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

    df['Charges'] = df.index.map(
        lambda peptide_id:
        charge_mapping.get(peptide_id, set()),
    )
    df['Masses'] = df.index.map(
        lambda peptide_id:
        mass_mapping.get(peptide_id, set()),
    )
    df['RTs'] = df.index.map(
        lambda peptide_id:
        rt_mapping.get(peptide_id, set()),
    )
    df['Intensities'] = df.index.map(
        lambda peptide_id:
        i_mapping.get(peptide_id, set()),
    )

    return df


def _set_defaults(df):
    df['Validated'] = False
    df['Fold Change'] = np.nan
    df['p-value'] = np.nan

    return df


def _get_filenames(df, cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        files = cursor.execute(
            '''
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
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        files = cursor.execute(
            '''
            SELECT
            TargetPsms.PeptideID,
            WorkflowInputFiles.PhysicalFileName

            FROM TargetPsms

            JOIN WorkflowInputFiles
            ON WorkflowInputFiles.FileID=TargetPsms.SpectrumFileID
            ''',
        )
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    mapping = {
        index: os.path.split(val)[-1]
        for index, val in files
    }

    df['Raw Paths'] = df.index.map(
        lambda peptide_id:
        mapping.get(peptide_id, {})
    )

    return df


def _get_q_values(df, cursor, pd_version):
    df['q-value'] = np.nan
    q_vals = []

    if pd_version[:2] in [(1, 4)]:
        fields = cursor.execute(
            '''
            SELECT
            CustomDataFields.FieldID,
            CustomDataFields.DisplayName

            FROM CustomDataFields
            ''',
        )

        field_ids = [
            field_id
            for field_id, name in fields
            if name in ['q-Value']
        ]

        if not field_ids:
            return df

        q_vals = cursor.execute(
            '''
            SELECT
            CustomDataPeptides.PeptideID,
            CustomDataPeptides.FieldValue

            FROM CustomDataPeptides

            WHERE CustomDataPeptides.FieldID IN ({})
            '''.format(
                ', '.join('?' * len(field_ids))
            ),
            field_ids,
        )
    elif pd_version[:2] in [(2, 2)]:
        try:
            q_vals = cursor.execute(
                '''
                SELECT
                TargetPsms.PeptideID,
                TargetPsms.PercolatorqValue

                FROM TargetPsms
                ''',
            )
        except sqlite3.OperationalError:
            pass
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    if q_vals:
        indices, vals = zip(*[
            (index, val)
            for index, val in q_vals
            if index in df.index
        ])

        df.at[indices, 'q-value'] = vals

    return df


def _is_pmod(mod):
    return (
        not mod.nterm and
        not mod.cterm and
        mod.mod_type in ['Phospho']
    )


def _sort_mods(mods):
    return tuple(
        sorted(
            mods,
            key=lambda x: (x.rel_pos, x.nterm, x.cterm, x.mod_type),
        )
    )


def _reassign_mods(mods, psp_val, probability_cutoff=75):
    reassigned = False
    ambiguous = False

    # phophoRS example format: 'T(4): 99.6; S(6): 0.4; S(10): 0.0'
    # Error messages include: 'Too many isoforms'
    if psp_val is None:
        psp_val = ''

    psp_val = [
        RE_PSP.match(i.strip())
        for i in psp_val.split(';')
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
    psp_val_f = [i for i in psp_val if i[2] > probability_cutoff]

    if len(p_mods) != len(psp_val_f):
        LOGGER.debug(
            'Not enough info to assign phophosite: {}'.format(psp_val)
        )
        ambiguous = True
    elif set(i.rel_pos + 1 for i in p_mods) != set(i[1] for i in psp_val_f):
        p_mods = [
            data_sets.Modification(
                rel_pos=i[1] - 1,
                mod_type='Phospho',
                nterm=False,
                cterm=False,
                sequence=p_mods[0].sequence,
            )
            for i in psp_val_f
        ]
        reassigned = True

        mods = data_sets.Modifications(
            mods=_sort_mods(o_mods + p_mods),
        )

        for mod in mods.mods:
            mod.sequence.modifications = mods

    return mods, reassigned, ambiguous


def _get_phosphors(df, cursor, pd_version, name=None):
    df['Ambiguous'] = False
    psp_vals = []

    if pd_version[:2] in [(1, 4)]:
        fields = cursor.execute(
            '''
            SELECT
            CustomDataFields.FieldID,
            CustomDataFields.DisplayName

            FROM CustomDataFields
            ''',
        )
        field_ids = [
            field_id
            for field_id, name in fields
            if name in ['phosphoRS Site Probabilities']
        ]

        if not field_ids:
            return df

        psp_vals = cursor.execute(
            '''
            SELECT
            CustomDataPeptides.PeptideID,
            CustomDataPeptides.FieldValue

            FROM CustomDataPeptides

            WHERE CustomDataPeptides.FieldID IN ({})
            '''.format(
                ', '.join('?' * len(field_ids))
            ),
            field_ids,
        )
    elif pd_version[:2] in [(2, 2)]:
        for ptmrs_col in [
            'ptmRSPhosphoSiteProbabilities',
            'ptmRSPhosphorylationSiteProbabilities',
        ]:
            try:
                psp_vals = cursor.execute(
                    '''
                    SELECT
                    TargetPsms.PeptideID,
                    TargetPsms.{}

                    FROM TargetPsms
                    '''.format(ptmrs_col),
                )
            except sqlite3.OperationalError:
                pass
            else:
                break
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    changed_peptides = 0

    for pep_id, psp_val in psp_vals:
        if pep_id not in df.index:
            continue

        old_mods = df.loc[pep_id]['Modifications']

        new_mods, reassigned, ambiguous = _reassign_mods(old_mods, psp_val)

        if reassigned:
            changed_peptides += 1

        df.at[pep_id, 'Modifications'] = new_mods
        df.at[pep_id, 'Ambiguous'] = ambiguous

    LOGGER.info(
        '{}: -- Reassigned {} phosphosites using phosphoRS'.format(
            name,
            changed_peptides,
        )
    )

    return df


def _get_species(cursor, pd_version):
    if pd_version[:2] in [(1, 4)]:
        fields = cursor.execute(
            '''
            SELECT
            ProcessingNodeParameters.ParameterValue

            FROM ProcessingNodeParameters

            WHERE ProcessingNodeParameters.ParameterName='Taxonomy'
            ''',
        )
    elif pd_version[:2] in [(2, 2)]:
        fields = cursor.execute(
            '''
            SELECT
            ProcessingNodeCustomData.CustomValue

            FROM ProcessingNodeCustomData

            WHERE ProcessingNodeCustomData.Name='FASTA database information'
            ''',
        )
        fields = [
            (line[len('Taxonomy:'):].strip(),)
            for i in fields
            for line in i[0].split('\n')
            if line.startswith('Taxonomy:')
        ]
    else:
        raise Exception(
            'Unsupported Proteome Discoverer Version: {}'.format(pd_version)
        )

    species = set()

    for name, in fields:
        # '. . . . . . . . . . . . . . . . Homo sapiens (human)'
        # => 'Homo sapiens'
        species.add(' '.join(name.strip('. ').split(' ')[:2]))

    return species


def read_discoverer_msf(basename, msf_path=None, pick_best_psm=False):
    '''
    Read a Proteome Discoverer .msf file.

    Converts file contents into a pandas DataFrame similar to what one would
    get by exporting peptides to .txt directly from Discoverer.

    Parameters
    ----------
    path : str
    pick_best_psm : bool, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    if msf_path is None:
        msf_path = os.path.join(
            paths.MS_SEARCHED_DIR,
            basename,
        )

    if not os.path.exists(msf_path):
        raise Exception('Search database does not exist: {}'.format(msf_path))

    name = os.path.splitext(basename)[0]

    LOGGER.info(
        '{}: Loading ProteomeDiscoverer peptides...'.format(name)
    )

    with sqlite3.connect(msf_path) as conn:
        cursor = conn.cursor()

        pd_version = _get_pd_version(cursor)

        _update_label_names(cursor, pd_version)

        # Get any N-terminal quantification tags
        tag_names = _get_quant_tags(cursor, pd_version)

        # Read the main peptide properties
        df = _read_peptides(conn, pd_version, pick_best_psm=pick_best_psm)

        LOGGER.info(
            '{}: -- Loading information for {} peptides'.format(
                name,
                df.shape[0],
            )
        )

        df = _get_proteins(df, cursor, pd_version)
        df = _extract_sequence(df)
        df = _extract_confidence(df)
        df = _extract_spectrum_file(df)
        df = _get_modifications(df, cursor, pd_version)
        df = _get_phosphors(df, cursor, pd_version, name=name)
        df = _get_q_values(df, cursor, pd_version)
        df = _get_ms_data(df, cursor, pd_version)
        df = _get_filenames(df, cursor, pd_version)
        df = _get_quantifications(df, cursor, pd_version, tag_names)

        df = _set_defaults(df)

        species = _get_species(cursor, pd_version)

    df['Scan Paths'] = basename

    df.reset_index(inplace=True, drop=True)

    LOGGER.info(
        '{}: Loaded {} peptides'
        .format(
            os.path.splitext(basename)[0],
            df.shape[0],
        )
    )

    return df, species

__all__ = [
    'read_discoverer_msf',
]