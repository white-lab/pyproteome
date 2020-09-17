'''
This module contains code for phosphorylation motif analysis.

It includes functions for discrete motif enrichment as well as generation of
motif logos. These logos can be generated locally (:func:`.logo.make_logo`) or
via automated hooks into online tools (:func:`.plogo.make_logo`,
:func:`.weblogo.make_logo`, :func:`.icelogo.make_logo`).
'''
from . import (
    icelogo, logo, motif, neighborhood, phosphosite, plogo, weblogo,
)

from .motif import (
    generate_n_mers,
)

__all__ = [
    'icelogo',
    'logo',
    'motif',
    'neighborhood',
    'phosphosite',
    'plogo',
    'weblogo',

    'generate_n_mers',
]
