# Change Log

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
