.. pyproteome documentation master file, created by
   sphinx-quickstart on Fri Jan 22 11:16:19 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyproteome's documentation!
======================================

pyproteome is a Python package for interacting with proteomics data.

It includes modules for loading, processing, and analyzing proteomics data
collected by mass spectometry. This functionality allows users to automatically
filter, normalize, and merge together data from proteome search files. It
analysis toolkit includes the ability to cluster peptides that show correlated
changes, and perform motif and pathway enrichment analysis to study cell
signaling events.

Currently it only supports analyzing ProteomeDiscoverer `.msf` search files.

This package is designed to be used within an interactive computational
environment, such as Jupyter Notebook, alongside other data analysis packages.
This allows scientists to create their analysis workflow within a reproducable
code environment.

Getting Started
===============

To start, you will need a Python environment. Python 2 and 3 are both supported.
Windows users may wish to use `Anaconda <https://www.anaconda.com/download/>`_
to manage their Python environment and provide pyproteome's dependencies.

  1. Install pyproteome::

      $ pip install pyproteome

  2. Open Python and load your data sets::

      from pyproteome import *
      # if using IPython:
      %import_all

      # with search data located in "Searched/"
      datas = data_sets.load_all_data(
          channels={
              "KO 1": "126",
              "KO 2": "127_N",
              "KO 3": "127_C",
              "Ctrl 1": "128_N",
              "Ctrl 2": "128_C",
              "Ctrl 3": "131",
          },
          groups={
              "Control": ["Ctrl 1", "Ctrl 2", "Ctrl 3"]},
              "KO": ["KO 1", "KO 2", "KO 3"],
          },
      )

  3. Analyze your data::

      # Show sequence motifs in upregulated set of peptides
      logo.make_logo(
          datas["CK-H1-pST"],
          {"asym_fold": 1.5},
      )

      # Cluster peptides into groups with similar changes
      data, y_pred = cluster.auto.auto_cluster(datas["CK-H1-pST"])

      # Show volcano plot of peptides enriched in cluster 1
      volcano.plot_volcano_filtered(
          datas["CK-H1-pST"],
          {"series": y_pred == 1},
      )

      # and show a table listing significantly changing peptides
      tables.changes_table(
          datas["CK-H1-pST"].filter(p=.001),
          sort="Fold Change",
      )

      # Perform Phospho Set Enrichment Analysis (PSEA)
      pathways.psea(
          datas["CK-H1-pST"],
          metric="zscore",
      )

      # Export the data set as an .csv file
      tables.write_csv(
        datas["CK-H1-pST"],
      )

Contents
========

.. toctree::
   :maxdepth: 5

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
