'''
This module provides functionality for accessing public data through
PRIDE PRoteomics IDEntifications (PRIDE) / Proteome Xchange.
'''

import os
# XXX: This should be a safer alternative package. Otherwise users could be
# vulnerable to a MITM attack
import xml.etree.ElementTree as ET

import requests


META_DATA_URL = (
    'http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={}&'
    'outputMode=XML&test=no'
)


def list_data_set(accession):
    '''
    Lists files contained in a deposition on PRIDE.

    Information is fetched from pride.META_DATA_URL.

    Parameters
    ----------
    accession : str

    Returns
    -------
    info_list : list of :class:`xml.etree.ElementTree`
        Information on files available in a repository.

    Examples
    --------
    >>> lst = pride.list_data_set('PXD003660')
    >>> lst[0].get('name')
    '20140524_MCF10A_E20VR1_ETP_TMT10.raw'
    >>> lst[0].find('cvParam').get('value')
    'ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD003660/20140524_MCF10A_E20VR1_ETP_TMT10.raw'
    '''
    assert accession.startswith('PXD')

    proj_id = int(accession[3:])

    # Fetch xml file using requests
    meta_data = requests.get(META_DATA_URL.format(proj_id))
    meta_data.raise_for_status()

    root = ET.fromstring(meta_data.text)

    return [
        i
        for i in root.findall('DatasetFileList/DatasetFile')
    ]


def fetch_data_set(accession, files=None):
    '''
    Fetches files from a deposition on PRIDE.

    Parameters
    ----------
    accession : str
        A PRIDE accession ID. i.e. 'PXD001038'
    files : dict of str, str
        Download individual files to a specific location. By default, this
        function downloads all files to the current working directory.

    Returns
    -------
    file_list : list of str
        Files downloaded from a repository.

    Examples
    --------
    >>> pride.fetch_data_set(
    ...     'PXD001038',
    ...     files={'HJ070512_OCTFF_B2_All5Fractions_PeptideSummary.zip': '.'},
    ... )
    ['HJ070512_OCTFF_B2_All5Fractions_PeptideSummary.zip']
    '''
    ds = list_data_set(accession)
    ret = []

    for file_root in ds:
        name = file_root.get('name')

        if files and name not in files:
            continue

        if isinstance(files, dict) and name in files:
            folder = files[name]
        else:
            folder = os.getcwd()

        out_path = os.path.join(folder, name)
        file_url = file_root.find('cvParam').get('value')

        # Requests cannot fetch FTP files
        if file_url.startswith('ftp://'):
            file_url = 'http://{}'.format(file_url[len('ftp://'):])

        response = requests.get(file_url, stream=True)
        response.raise_for_status()

        with open(out_path, mode='wb') as f:
            for block in response.iter_content(1024):
                f.write(block)

        ret.append(name)

    return ret
