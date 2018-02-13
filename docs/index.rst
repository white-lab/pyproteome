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

  1. Install pyproteome: ``pip install pyproteome``

  2. Load your data sets:

  ``
  from pyproteome import *
  # if using IPython:
  %import_all

  # with search data located in "Searched/"
  datas = data_sets.load_all_data(
      channels={
          "126": "3130 CK",
          "127": "3131 CK-p25",
          "128": "3145 CK-p25",
          "129": "3136 CK-p25",
          "130": "3148 CK",
          "131": "3157 CK",
      },
      groups=[
          {"CK-p25": ["3131 CK-p25", "3145 CK-p25", "3146 CK-p25"]},
          {"CK": ["3130 CK", "3148 CK", "3157 CK"]},
      ],
  )
  ``

  3. Analyze your data:
  
  ``
  # Show sequence motifs in upregulated set of peptides
  logo.make_logo(
      datas["CK-H1-pST"],
      {"asym_fold": 1.5, "p": 0.01},
  )

  # Cluster peptides into groups with similar changes
  data, y_pred = cluster.auto.auto_cluster(datas["CK-H1-pST"])

  # Show volcano and sequence motif logo for peptides enriched in cluster 1
  volcano.plot_volcano_filtered(
      datas["CK-H1-pST"],
      {"series": y_pred == 1},
  )
  logo.make_logo(
      datas["CK-H1-pST"],
      {"series": y_pred == 1},
  )
  ``

Python 2 and 3 are both supported. Windows users may wish to use Anaconda to
manage their Python environment and provide pyproteome's dependencies.

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

