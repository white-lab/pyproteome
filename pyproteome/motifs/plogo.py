# -*- coding: utf-8 -*-
from __future__ import division

# Built-ins
import logging
import math
import requests
import time

try:
    from IPython.display import Image
except ImportError:
    Image = None

from . import motif


LOGGER = logging.getLogger('pyproteome.plogo')
PLOGO_BASE = 'https://plogo.uconn.edu/main'


def format_title(
    f=None,
    data=None,
):
    '''
    Generates a title automatically from a given data set and list of filters.

    Parameters
    ----------
    f : dict or list of dict
    data : :class:`pyproteome.data_sets.data_set.DataSet`

    Returns
    -------
    str
    '''
    title = []

    if 'fold' in f:
        title.append(
            'abs(FC) > {:.2f}'.format(f['fold'])
        )

    if 'asym_fold' in f:
        title.append(
            'FC {} {:.2f}'.format(
                '>' if f['asym_fold'] > 1 else '<',
                f['asym_fold'],
            )
        )

    if 'p' in f:
        title.append(
            'p < {:.2f}'.format(f['p'])
        )

    title = ', '.join(title)

    if data:
        title = '{} - {}'.format(data.name, title)

    return title


def make_logo(data, f, **kwargs):
    '''
    Wraps :func:`.plogo`, generating the list of foreground and background
    peptide sequences from a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.data_set.DataSet`
    f : dict or list of dict
        Argument passed to :func:`pyproteome.data_sets.data_set.DataSet.filter`.
    kwargs : dict
        Arguments passed to :func:`.plogo`.

    Returns
    -------
    str or :class:`IPython.display.Image`
    '''
    nmer_args = motif.get_nmer_args(kwargs)

    fore = [
        n
        for n in motif.generate_n_mers(
            data.filter(f)['Sequence'],
            **nmer_args
        )
    ]

    back = [
        n
        for n in motif.generate_n_mers(
            data['Sequence'],
            **nmer_args
        )
    ]
    return plogo(
        fore, back,
        title=format_title(data=data, f=f),
        **kwargs
    )


def _plogo_wait_job(s, job, delay=1):
    while True:
        response = s.get(
            '{}/getPlogo/{}'.format(PLOGO_BASE, job),
        )
        response.raise_for_status()
        json = response.json()

        if json['complete']:
            break

        time.sleep(delay)


def _check_plogo_response(response, message=''):
    response.raise_for_status()
    json = response.json()
    assert json['message'] == message
    assert json['success']
    return json


def plogo(
    foreground, background,
    fix_letter_pos=None,
    title='',
    width=800,
    height=600,
    ymax=None,
):
    '''
    Wraps calls to the pLogo web server [1]_, returning an image showing the enrichment
    of a sequence in a foreground set compared to a background set.

    Parameters
    ----------
    foreground : list of str
    background : list of str
    fix_letter_pos : list of tuple of (str, int), optional
    title : str, optional
    width : int, optional
    height : int, optional
    ymax : float, optional

    Returns
    -------
    str or :class:`IPython.display.Image`

    Notes
    -----
    .. [1] O’Shea, Joseph P et al. “pLogo: A Probabilistic Approach to
       Visualizing Sequence Motifs.” Nature Methods 10.12 (2013): 1211–1212.
       Web.
    '''
    if fix_letter_pos is None:
        fix_letter_pos = []

    assert len(foreground) > 0
    assert len(background) > 0

    letter_width = len(background[0])

    s = requests.Session()

    response = s.post(
        '{}/uploadforeground/'.format(PLOGO_BASE),
        data={
            'fg': '\n'.join(foreground),
            'fix': 1,
            'dataType': 'text',
            'remove_duplicates': True,
        },
    )
    json = _check_plogo_response(
        response, message='Foreground uploaded. See stats.',
    )
    fore_dir = json['dir']

    response = s.post(
        '{}/uploadbackground/'.format(PLOGO_BASE),
        data={
            'bg_text': '\n'.join(background),
            'dataType': 'text',
            'datasetType': None,
            'name': 'unknown',
            'remove_duplicates': True,
        },
    )
    json = _check_plogo_response(response)
    back_dir = json['dir']
    assert fore_dir == back_dir

    response = s.post(
        '{}/initplogo/'.format(PLOGO_BASE),
        data={
            'width': letter_width,
            'fixedLetter': 'S',
            'fixedPosition': 7,
            'fix': 0,
            'jobName': 'pyproteome.motif.plogo',
            'subtract_fg': False,
            'remove_duplicates': True,
            'from_ksdb': False,
        },
    )

    json = _check_plogo_response(response)
    job = json['state']['jobId']
    _plogo_wait_job(s, job)

    for letter, pos in fix_letter_pos:
        response = s.post(
            '{}/initplogo/'.format(PLOGO_BASE),
            data={
                'width': letter_width,
                'fixedLetter': letter,
                'fixedPosition': math.floor(letter_width / 2) + pos,
                'fix': 1,
                'jobName': 'pyproteome.motif.plogo - fix',
                'subtract_fg': False,
                'remove_duplicates': True,
                'from_ksdb': False,
                'referenceId': job,
            },
        )

        json = _check_plogo_response(response)
        job = json['state']['jobId']
        _plogo_wait_job(s, job)

    response = s.post(
        '{}/redrawplogo/'.format(PLOGO_BASE),
        data={
            'options[title]': title,
            'options[width]': width,
            'options[height]': height,
            'options[plogoMode]': 'absolute',
            'options[maxAbsoluteValue]': ymax,
            'jobId': job,
        },
    )
    json = _check_plogo_response(response)
    job = json['state']['jobId']
    _plogo_wait_job(s, job)

    png = json['state']['job']['image_url']

    if Image:
        return Image(url=png)
    else:
        return png
