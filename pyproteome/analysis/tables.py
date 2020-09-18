
from __future__ import absolute_import, division

import collections
import os
import re

import pandas as pd

import pyproteome as pyp


def _prep_csv(data=None, postfix='table', folder_name=None, csv_name=None):
    if csv_name is None:
        csv_name = '{}.csv'.format(
            postfix,
        )

    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub='Tables',
    )

    return os.path.join(folder_name, csv_name)


def _get_table_title(f=None, running_title=None):
    if running_title is None:
        running_title = []

    if f is not None:
        if 'asym_fold' in f:
            running_title.append(
                'Upregulated' if f['asym_fold'] > 1 else 'Downregulated'
            )

        if 'p' in f:
            running_title.append(
                'p-{:.3e}'.format(f['p'])
            )

        if 'group_a' in f or 'group_b' in f:
            running_title.append(
                '{}vs{}'.format(
                    f.get('group_a', ''),
                    f.get('group_b', ''),
                )
            )

    return '-'.join(running_title)


def motif_table(
    data, f,
    p=0.05,
    sort='p-value',
    **kwargs
):
    '''
    Run a motif enrichment algorithm on a data set and display the
    significantly enriched motifs.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    f : dict or list of dict
    p : float, optional
    sort : str, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`

    See Also
    --------
    :func:`pyproteome.motifs.motif.run_motif_enrichment`
    '''
    hits = pyp.motifs.motif.run_motif_enrichment(
        data, f,
        **kwargs
    )[0]

    hits = hits.sort_values(
        sort,
        ascending=True if sort == 'p-value' else False,
    )

    hits = hits[
        hits[
            'pp-value' if kwargs.get('pp_value', False) else 'p-value'
        ] < p
    ]

    return hits.style.set_table_styles([
        {'selector': '*', 'props': [('font-family', 'monospace')]},
        {'selector': 'th:first-child', 'props': [('display', 'none')]},
    ])


def changes_table(
    data,
    sort='p-value',
):
    '''
    Show a table of fold changes and p-values for each unique peptide in a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    sort : str, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    psms = getattr(data, 'psms', data)

    psms = psms[
        [
            'Proteins', 'Sequence', 'Modifications',
            'Fold Change', 'p-value', 'Validated',
        ]
    ].copy()
    psms['Sequence'] = psms['Sequence'].apply(
        lambda x: '{} ({})'.format(x, x.modifications)
    )
    psms['Uniprot Accessions'] = psms['Proteins'].apply(
        lambda x: '; '.join(x.accessions)
    )

    if sort == 'Fold Change':
        psms['Fold Change-Sort'] = psms['Fold Change'].apply(
            lambda x: max([x, 1 / x])
        )
        psms.sort_values('Fold Change-Sort', inplace=True, ascending=False)
        psms.drop('Fold Change-Sort', axis=1, inplace=True)
    else:
        psms.sort_values(sort, inplace=True, ascending=True)

    psms.drop('Modifications', axis=1, inplace=True)

    # back_colors = {
    #     True: '#BBFFBB',  # light green
    #     False: '#FFBBBB',  # light red
    # }

    if psms.empty:
        return psms

    # return psms.style.apply(  # Color validated rows
    #     lambda row: [
    #         'background-color: ' + back_colors[row['Validated']]
    #         for _ in row
    #     ],
    #     axis=1,
    # )
    return psms.style.set_table_styles(  # Hide index and 'Validated' columns
        [
            {'selector': 'th:first-child', 'props': [('display', 'none')]},
            {'selector': 'td:last-child', 'props': [('display', 'none')]},
            {'selector': 'th:last-child', 'props': [('display', 'none')]},
            {'selector': '*', 'props': [('text-align', 'left')]},
        ]
    )


def ptmsigdb_changes_table(
    data,
    sort='p-value',
    folder_name=None,
    csv_name=None,
):
    '''
    Show a table of fold changes and p-values for PTMSigDB.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    sort : str, optional
    folder_name : str, optional
    csv_name : str, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    csv_name = _prep_csv(
        data=data,
        folder_name=folder_name,
        csv_name=csv_name,
        postfix=_get_table_title(running_title=['ptmsigdb']),
    )

    psms = getattr(data, 'psms', data).copy()
    psms = psms.dropna(subset=('Fold Change',))

    psms['Protein Description'] = psms['Proteins'].apply(
        lambda x: x.proteins[0].description
    )
    psms['Gene'] = psms['Proteins'].apply(
        lambda x: x.genes[0]
    )
    psms['Uniprot Accession'] = psms['Proteins'].apply(
        lambda x: x.accessions[0]
    )
    psms['All Modifications'] = psms['Modifications']
    psms['Phospho Modifications'] = psms['Modifications'].apply(
        lambda x: x.get_mods([(None, 'Phospho')]).__str__(prot_index=0)
    )

    psms = psms[
        [
            'Protein Description',
            'Gene',
            'Uniprot Accession',
            'Sequence',
            'All Modifications',
            'Phospho Modifications',
            'Fold Change',
        ]
    ].copy()
    psms.sort_values('Fold Change', inplace=True, ascending=False)

    if csv_name:
        psms.to_csv(
            csv_name,
            index=False,
        )

    if psms.empty:
        return psms

    return psms.style.set_table_styles(  # Hide index and 'Validated' columns
        [
            {'selector': 'th:first-child', 'props': [('display', 'none')]},
            {'selector': '*', 'props': [('text-align', 'left')]},
        ]
    )


def _ds_to_df(data, save_cols=None, sample_values=True):
    channels = [
        (name, data.channels[name])
        for name in data.samples
    ] if sample_values else []
    
    if save_cols is None:
        save_cols = ['Fold Change', 'p-value']
    else:
        save_cols = [i for i in save_cols if i in data.psms.columns]
    
    cols = [
        'Proteins', 'Sequence', #'Scan',
    ] + save_cols + [
        chan
        for _, chan in channels
    ]
    df = data.psms[cols].copy()

    df.rename(
        columns={
            chan: name
            for name, chan in channels
        },
        inplace=True,
    )
    df.insert(
        1, 'Genes',
        df['Proteins'].apply(
            lambda x:
            ' / '.join(x.genes)
        ),
    )
    df.insert(
        2, 'Uniprot Accessions',
        df['Proteins'].apply(
            lambda x:
            ' / '.join(x.accessions)
        ),
    )
    df.insert(
        3, 'Modifications',
        df['Sequence'].apply(
            lambda x: str(x.modifications)
        ),
    )
    df['Sequence'] = df['Sequence'].apply(str)
    # df['Scan'] = df['Scan'].apply(
    #     lambda x:
    #     ', '.join([str(i) for i in x])
    #     if isinstance(x, collections.Iterable) else
    #     str(x)
    # )

    if 'p-value' in cols:
        df.sort_values('p-value', inplace=True, ascending=True)

    return df


def write_csv(data, folder_name=None, out_name='DataSet.csv'):
    '''
    Write information for a single data set to a .csv file.

    Sheets are populated with protein, peptide, scan, and quantification values
    for all peptide-spectrum matches contained within a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    folder_name : str, optional
    out_name : str, optional

    Returns
    -------
    path : str
        Path to .xlsx file.
    '''

    out_name = _prep_csv(
        data=None,
        folder_name=folder_name,
        csv_name=out_name,
    )

    df = _ds_to_df(data)

    df.to_csv(
        out_name,
        index=False,
    )

    return out_name


def write_full_tables(
    datas, 
    save_cols=None, 
    sample_values=True,
    folder_name=None, 
    out_name='Full Data.xlsx', 
):
    '''
    Write information for a list of data sets to sheets of a .xlsx file.

    Sheets are populated with protein, peptide, scan, and quantification values
    for all peptide-spectrum matches contained within a data set.

    Parameters
    ----------
    datas : list of :class:`pyproteome.data_sets.DataSet`
    save_cols : list of str, optional
        Extra column names to save from in each dataset.
    sample_values : bool, optional
        Save normalized TMT values for each sample to the output.
    folder_name : str, optional
    out_name : str, optional

    Returns
    -------
    path : str
        Path to .xlsx file.
    '''

    out_name = _prep_csv(
        data=None,
        folder_name=folder_name,
        csv_name=out_name,
    )

    writer = pd.ExcelWriter(out_name, engine='xlsxwriter')

    for data in datas:
        df = _ds_to_df(
            data, 
            save_cols=save_cols,
            sample_values=sample_values,
        )

        ws_name = re.sub(
            '/',
            '+',
            data.name,
        )
        df.to_excel(
            writer,
            sheet_name=ws_name,
            index=False,
        )

        ws = writer.sheets[ws_name]
        
        ws.freeze_panes(1, 0)
        
        ws.set_column(0, 0, 60)
        ws.set_column(1, 1, 15)
        ws.set_column(2, 2, 15)
        ws.set_column(3, 3, 20)
        ws.set_column(4, 4, 20)

        ws.autofilter(0, 0, df.shape[0], df.shape[1] - 1)

        for ind, col in enumerate(df.columns):
            args = 0, ind, df.shape[0], ind

            if col.endswith('-p'):
                args += {'type': '2_color_scale', 'minimum': 0, 'maximum': 5e-2},
            elif col.endswith('-FC'):
                args += {'type': '3_color_scale'},
            elif col.endswith('-Corr') or col.endswith('-Correlation'):
                args += {'type': '3_color_scale', 'minimum': -1, 'maximum': 1},
            else:
                continue

            ws.conditional_format(*args)

    writer.save()

    return out_name
