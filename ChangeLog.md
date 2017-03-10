# Change Log

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
