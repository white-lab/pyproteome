from __future__ import division

# Built-ins
import logging
import requests

try:
    from IPython.display import Image
except ImportError:
    Image = None

from . import motif


LOGGER = logging.getLogger('pyproteome.icelogo')
ICELOGO_BASE = 'http://iomics.ugent.be/icelogoserver'


def make_logo(data, f, m=None, letter_mod_types=None, **kwargs):
    nmer_args = motif.get_nmer_args(kwargs)

    fore = [
        n
        for n in motif.generate_n_mers(
            data.filter(f)['Sequence'],
            **nmer_args
        )
        if not m or m.match(n)
    ]

    back = [
        n
        for n in motif.generate_n_mers(
            data['Sequence'],
            **nmer_args
        )
    ]
    return icelogo(
        fore, back,
        **kwargs
    )


def icelogo(
    foreground, background,
    title='',
    width=800,
    height=600,
    pvalue=0.05,
    scoring='foldChange',
):
    '''
    Wraps calls to iceLogo [1]_, returning an image showing the enrichment of a
    sequence in a foreground set compared to a background set.

    Parameters
    ----------
    foreground : list of str
    background : list of str
    title : str, optional
    width : int, optional
    height : int, optional
    pval : float, optional
    scoring : string, optional

    Returns
    -------
    fig : :class:`IPython.display.Image`

    Notes
    -----
    .. [1] Colaert, N., Helsens, K., Martens, L., Vandekerckhove, J., &
        Gevaert, K. (2009). Improved visualization of protein consensus
        sequences by iceLogo. Nature Methods, 6(11), 786–787.
        http://doi.org/10.1038/nmeth1109-786
    '''
    if len(foreground) < 1 or len(background) < 1:
        return None

    assert scoring in ('percentage', 'foldChange')

    letter_width = len(background[0])

    s = requests.Session()

    response = s.post(
        '{}/data/logo'.format(ICELOGO_BASE),
        data={
            'positiveSequences': '\n'.join(foreground),
            'reference': 'reference_set',
            'negativeSequences': '\n'.join(background),
            'start': -(letter_width // 2),
            'visualisationType': 'iceLogo',
            'species': '',
            'colors[]': [
                '2b00e3', 'da00ba', 'd20012', 'bcc800', 'cd6300', '0000e5',
                'c9af00', 'd60064', 'cf3d00', 'd8008f', 'cb8900', 'd4003b',
                '8000df', 'd11500', '94c600', '5600e1', 'd200db', 'a900dd',
                '47c200', '21c000', '6dc400',
            ],
            'width': width,
            'height': height,
            'pValue': pvalue,
            'scoringSystem': scoring,
            'aaMatrix': '',
            'substitutionMatrix': '',
        },
    )
    response.raise_for_status()
    path = response.content.decode('utf-8')

    png = '{}/generated_images/{}.png'.format(ICELOGO_BASE, path)

    if Image:
        return Image(url=png)
    else:
        return png
