'''
This module provides functionality for manipulating proteomics data sets.

Functionality includes merging data sets and interfacing with attributes in a
structured format.
'''

# Built-ins
from __future__ import absolute_import, division

from pyproteome import levels

import logging

# Core data analysis libraries
import numpy as np


LOGGER = logging.getLogger('pyproteome.constand')
CONSTAND_METHODS = {
    'mean': np.nanmean,
    'median': np.nanmedian,

    # Fit a gaussian function to each row/column and select its max point
    'kde': lambda k, **kw: np.apply_along_axis(levels.kde_max, kw['axis'], k),
}
'''
Methods to estimate the center of a data set's rows or columns.

One of ['mean', 'median', 'kde']. 'mean': applies :func:`numpy.nanmean`,
'median' applies :func:`numpy.nanmedian`, and 'kde' applies
:func:`.levels.kde_max`.
'''

CONSTAND_ERR_METHODS = {
    'mean': lambda x, **kw: np.sqrt(abs(np.nanmean(x, **kw) - 1)).sum() / 2,
    'median': lambda x, **kw: np.sqrt(abs(np.nanmedian(x, **kw) - 1)).sum() / 2,
    'kde': lambda x, **kw: np.sqrt(abs(np.apply_along_axis(levels.kde_max, kw['axis'], x) - 1)).sum() / 2,
}
DEFAULT_CONSTAND_ROW = 'mean'
DEFAULT_CONSTAND_COL = 'median'


def constand(
    ds,
    inplace=False,
    n_iters=25, 
    tol=1e-5,
    row_method=None,
    col_method=None,
):
    '''
    Normalize channels for intra-run comparisons. Iteratively fits the matrix
    of quantification values such that each row and column are centered around
    a calculated value. Uses row means and column median values for centering
    by default. See `constand.CONSTAND_METHODS` for other options.

    Parameters
    ----------
    ds : :class:`.data_sets.DataSet`
        Data set to apply CONSTANd normalization on.
    inplace : bool, optional
        Modify this data set in place.
    n_iters : int, optional
        Max number of normalization iterations. Rows are normalized on
        the odd step and columns are normalized on the even step.
    tol : float, optional
        Minimum error tolerance to use to end iterations early.
    row_method : str, optional
        Row normalization method to use. Default value is 'mean'.
    col_method : str, optional
        Column normalization method to use. Default value is 'median'.

    Returns
    -------
    ds : :class:`.data_sets.DataSet`
    '''
    new = ds

    if row_method is None:
        row_method = DEFAULT_CONSTAND_ROW

    if col_method is None:
        col_method = DEFAULT_CONSTAND_COL

    if not inplace:
        new = new.copy()

    channels = list(new.channels.values())

    new_channels = new.channels

    k = new[channels].values
    err = np.inf
    
    row_fn = lambda x: CONSTAND_METHODS[row_method](x, axis=1)
    col_fn = lambda x: CONSTAND_METHODS[col_method](x, axis=0)

    row_err_fn = lambda x: CONSTAND_ERR_METHODS[col_method](x, axis=0)
    col_err_fn = lambda x: CONSTAND_ERR_METHODS[row_method](x, axis=1)

    for ind in range(1, n_iters + 1):
        if ind % 2 == 1:
            # In the odd step the rows are fitted to match the row marginals
            # (i.e. constraints):
            #
            # K^(2t + 1) = R^(t + 1) * K ^ (2t)
            #
            # The row multipliers R^(t + 1) are computed such that the mean
            # of the reporter ion intensities equals 1:
            r = 1 / row_fn(k)

            # k^(2t + 1)
            k = np.einsum('..., ...', r, k.T).T
            err = row_err_fn(k)
        else:
            # In the even step the columns are fitted to match the column
            # marginals (i.e. constraints):
            #
            # K^(2t + 2) = K^(2t + 1) * S^(t + 1)
            #
            # The column multipliers S^(t + 1) are computed such that the
            # mean of the reporter ion intensities equals 1:
            s = 1 / col_fn(k)

            # k^(2t + 2)
            k = np.einsum('..., ...', k, s)
            err = col_err_fn(k)

        # print(ind, err)

        if err < tol:
            break

    LOGGER.info(
        '{}Applied CONSTANd normalization: iters = {}, err = {:.2e}.'.format(
            new.name + ': ' if new.name else '',
            ind,
            err,
        )
    )

    for i, (key, norm_key) in enumerate(zip(
        new.channels.values(),
        new_channels.values(),
    )):
        new.psms[norm_key] = k[:, i]

        if key != norm_key:
            del new.psms[key]

    new.channels = new_channels
    new.groups = new.groups.copy()

    new.update_group_changes()

    return new
