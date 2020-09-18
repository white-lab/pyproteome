import os
import math
import requests

try:
    from IPython.display import Image
except ImportError:
    Image = None

from . import motif


WEBLOGO_URL = 'http://weblogo.threeplusone.com/create.cgi'


def make_logo(
    data,
    **kwargs
):
    '''
    Create a sequence logo figure.

    Logos are created based on the frequencies of peptides in a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    '''

    nmer_args = motif.get_nmer_args(kwargs)

    nmers = motif.generate_n_mers(
        data.psms['Sequence'],
        **nmer_args
    )
    nmer_txt = '\n'.join(nmers) + '\n'
    first_index = math.ceil(-nmer_args.get('n', 15) / 2)

    response = requests.post(
        WEBLOGO_URL,
        files=(
            ('sequences_file', ('', nmer_txt, 'application/octet-stream')),
            ('sequences_url', (None, '')),
            ('sequences', (None, '')),
            ('cmd_create', (None, '')),
            ('logo_title', (None, data.name)),
            ('format', (None, 'png_print')),
            ('alphabet', (None, 'alphabet_auto')),
            ('stack_width', (None, 'large')),
            ('stacks_per_line', (None, '40')),
            ('unit_name', (None, 'bits')),
            ('first_index', (None, str(first_index))),
            ('logo_start', (None, '')),
            ('logo_end', (None, '')),
            ('scale_width', (None, 'true')),
            ('composition', (None, 'comp_auto')),
            ('percentCG', (None, '')),
            ('show_errorbars', (None, 'true')),
            # ('show_fineprint', (None, 'true')),
            ('show_xaxis', (None, 'true')),
            ('xaxis_label', (None, '')),
            ('show_yaxis', (None, 'true')),
            ('yaxis_label', (None, '')),
            ('yaxis_scale', (None, 'auto')),
            ('yaxis_tic_interval', (None, '1.0')),
            ('color_scheme', (None, 'color_auto')),
            ('symbols0', (None, '')),
            ('symbols1', (None, '')),
            ('symbols2', (None, '')),
            ('symbols3', (None, '')),
            ('symbols4', (None, '')),
            ('color0', (None, '')),
            ('color1', (None, '')),
            ('color2', (None, '')),
            ('color3', (None, '')),
            ('color4', (None, '')),
        ),
        stream=True,
    )

    response.raise_for_status()

    if Image:
        return Image(response.content, format='png')
    else:
        return response.content
