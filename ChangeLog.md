# Change Log

### 0.9.1 (2019-01-10)

Bug fixes

  - Fixed minor bugs in figure layouts.

### 0.9.0 (2019-01-10)

Features

  - Added support for ProteomeDiscoverer 2.2.
  - Improved the layout of several analysis figures.
  - Added cell type, gene ontology enrichment functions.

## 0.8.0 (2018-08-27)

Changes

  - Reorganized DataSet, Protein, Sequence, Modification modules into a single
    package (pyproteome.data_sets).

## 0.7.2 (2018-08-25)

Features

  - Significant improvements to documentation.
  - Cleaned up several figure layouts.

Bug fixes

  - Fixed install from pip.

## 0.7.1 (2018-08-22)

Features

  - More consistent figure + text sizes.
  - Added tables.write_csv function.

## 0.7.0 (2018-08-20)

Features

  - Many improvements to GSEA / PSEA calculations. Added functionality for
    remapping phosphosites between species.
  - Renamed several functions / parameters for more consistent naming.
  - Added QC module for validation of peptide quantification robustness.
  - Re-wrote several plotting functions to use seaborn for pretty figures.
  - Changed normalization calculations to fit to the mode, rather than the
    median channel ratio.

## 0.6.1 (2018-03-22)

Bug fixes

  - Fixed PyPI deployment issues.

## 0.6.0 (2018-03-22)

Features

  - Improvements to ProteomeDiscoverer feature support (percolator, phosphoRS).
  - Many organizational changes.
  - Improvements to clustering module.
  - Added pathways module for GSEA / PSEA.

## 0.5.0 (2018-01-20)

Features

  - Completed clustering module.
  - Added much more extensive integration testing, fixing many py2 bugs.
  - Re-organized analysis into package.
  - Added code for normalizing across data sets with no shared channels.
  - Cleared out old code for reading tab-delimited Discoverer exports.
  - Save plots into organized folders.

## 0.4.0 (2018-01-13)

Features

  - Added pick_best_ptm option to DataSet, enabled by default. This selects
    only the peptide with the best PTM assignment to be used for each scan. It
    is recommended to disable this feature on large pS/T data sets.
  - Sequence.\_\_contains\_\_ supports string arguments now.
  - Added individual sample quantitations to write_full_tables()
  - Automatically apply inter-normalization when merging multiple data sets.
  - Filter by several new parameters.
  - Added motifs, cluster, brainrnaseq modules.
  - Added integration tests.
  - Reorganized several analysis functions.
  - Removed unused / untested enrichment module.

Bugfixes

  - Fixed several runtime warnings from matplotlib / numpy / pandas.
  - Fixed many py2 issues.

## 0.3.2 (2017-03-10)

Bugfixes

  - Fixed levels.get_channel_levels not displaying all histogram plots
  - Cleared up clutter on analysis.volcano_plot

## 0.3.1 (2017-02-28)

Features

  - Drop notebook requirement down to just ipython.

Bugfixes

  - Fixed truncation of x/y-axis ticks on volcano plots.

## 0.3.0 (2017-02-28)

Features

  - Added support for > 2 groups.
  - Added support for inter-run normalization.

Bugfixes

  - Too many bug-fixes to count.

## 0.2.3 (2016-09-12)

Bugfixes

  - Fixed uniprot needing to be installed separately.

## 0.2.2 (2016-09-12)

Features

  - Added in uncommitted BCA assay changes.

## 0.2.1 (2016-08-18)

Features

  - Added CLI front end to pycamv (`python -m pycamv`)
  - Added pycamv support for isolation windows.

## 0.2.0 (2016-08-10)

Features

  - Filled in pycamv package, can now export full .camv files for
    [CAMV](https://github.com/naderm/CAMV) validation.
  - Added pyprotome.pride module for fetching public data sets.
  - Added pyproteome.discoverer module for reading Proteome Discoverer files.
  - Color pyproteome.analysis.snr_table() by peptide validation.
  - Show "Validated" column in pyproteome.analysis.write_full_tables().

Bugfixes

  - Fixed several bugs in calculating fragment ion masses.
  - Install all packages from pypi.

## 0.1.0 (2016-03-29)

Features

  - Added pycamv package, provides functionality for validating data sets.
  - Added pyproteome.bca module, provides functionality for processing BCA
    protein concentration measurements.
  - Added several new functions to pyproteome.analysis.
  - Added pyproteome.paths module for custom directory structures.

Bugfixes:

  - Fixed numerous py2 vs. py3 bugs.
  - Updated many module and function docstrings.
  - Added doc interlinks for function parameter types.

## 0.0.1 (2016-03-03)

Initial release
