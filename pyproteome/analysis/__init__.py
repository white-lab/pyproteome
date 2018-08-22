# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence
levels.
"""

from . import correlation, heatmap, plot, tables, qc, volcano

__all__ = [
    "correlation",
    "heatmap",
    "plot",
    "tables",
    "qc",
    "volcano",
]
