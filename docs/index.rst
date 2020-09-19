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

To start, you will need a Python environment. Python version >= 3.6 is recommended.
Windows users may wish to use `Anaconda <https://www.anaconda.com/download/>`_
to manage their Python environment and provide pyproteome's dependencies.

  1. Install pyproteome::

      $ pip install pyproteome

  2. Open Python and load your data sets, using :func:`pyproteome.data_sets.load_all_data`::

      from pyproteome import *
      # if using IPython:
      %import_all

      # Define sample:TMT channel mapping for each TMT-10plex analysis
      run1_chans = {
          'KO 1': '126',
          'KO 2': '127N',
          'KO 3': '127C',
          'KO 4': '128C',
          'KO 5': '128C',
          'Ctrl 1': '129N',
          'Ctrl 2': '129C',
          'Ctrl 3': '130N',
          'Ctrl 4': '130C',
          'Pooled Samples': 'Pooled Samples',
      }
      run2_chans = {
          'KO 6': '126',
          'KO 7': '127N',
          'KO 8': '127C',
          'KO 9': '128C',
          'KO 10': '128C',
          'Ctrl 5': '129N',
          'Ctrl 6': '129C',
          'Ctrl 7': '130N',
          'Ctrl 8': '130C',
          'Pooled Samples': 'Pooled Samples',
      }

      # This example demostrates loading and processing 6 different proteomics analyses
      # with 2 different TMT-10plex labeled samples, 3 different enrichment methods, and
      # 1 common pooled sample.
      #
      # ProteomeDiscoverer .msf search files are first stored in 'Searched/' as:
      #     Run1_pY.msf
      #     Run1_pSQTQ.msf
      #     Run1_Sup.msf
      #     Run2_pY.msf
      #     Run2_pSQTQ.msf
      #     Run2_Sup.msf
      datas = data_sets.load_all_data(
          # Assign run1 channels to all search files beginning with 'Run1_'
          # and run2 channels to all search files beginning with 'Run2_'
          chan_mapping={
              'Run1_': run1_chans,
              'Run2_': run2_chans,
          },
          
          # Apply CONSTANd normalization
          norm_mapping='constand',

          # Filter out peptides that do not pass the quality-control cutoffs
          filter_bad={
              'ion_score': 15,
              'isolation': 30,
              'median_quant': 1.5e3,
              'q': 1e-2,
          }

          # Merge pY, pSQ/pTQ, and global supernatant runs together,
          # then merge Run1 and Run2 together, normalized against their
          # common channel
          merge_mapping={
              'Run1': ['Run1_pY', 'Run1_pSQTQ', 'Run1_Sup'],
              'Run2': ['Run2_pY', 'Run2_pSQTQ', 'Run2_Sup'],
              'Merged': ['Run1', 'Run2],
          }

          # Create a list for each comparison group
          groups={
              'KO': [
                  'KO 1', 'KO 2', 'KO 3', 'KO 4', 'KO 5',
                  'KO 6', 'KO 7', 'KO 8', 'KO 9', 'KO 10',
              ],
              'Control': [
                  'Ctrl 1', 'Ctrl 2', 'Ctrl 3', 'Ctrl 4',
                  'Ctrl 5', 'Ctrl 6', 'Ctrl 7', 'Ctrl 8',
              ],
              'Pooled': [
                  'Pooled Samples',
              ],
          },
      )

  3. Analyze your data, using :func:`changes_table<pyproteome.analysis.tables.changes_table>`,
     :func:`make_logo<pyproteome.motifs.logo.make_logo>`,
     :func:`plot_volcano<pyproteome.analysis.volcano.plot_volcano>`,
     and :func:`psea<pyproteome.pathways.psea>`::

      # Show a table listing significantly changing peptides
      display(
          tables.changes_table(
              # fold -> Average Group Fold Change (FC): FC > 1.5 or FC < 1/1.5
              # p -> 2-sample t-test p-value between groups: p < 1e-2
              datas['Merged'].filter(fold=1.25, p=1e-2),

              # Sort by fold change, otherwise sort by p-value by default
              sort='Fold Change',
          )
      )

      # Show phosphorylation motifs in upregulated set of peptides
      logo.make_logo(
          datas['Merged'],
          {'asym_fold': 1.5, 'p': 1e-2},
      )

      # Show volcano plot of peptides enriched in cluster 1
      volcano.plot_volcano(
          datas['Merged'],
          fold=1.5,
          p=1e-3,
      )

      # Perform Phospho Set Enrichment Analysis (PSEA)
      pathways.psea(
          datas['Merged'],
          min_hits=15,
          pval=True,
          metric='zscore',
          p=.75,
          p_iter=500,
          max_pval=1e-2,
          max_qval=.25,
          n_cpus=4,
      )

      # Export the data set with fold changes as a .csv file
      tables.write_csv(
        datas['Merged'],
        out_name='Merged.csv',
      )

      # Export all quantification data from each data set to an excel table
      tables.write_full_tables(
        datas,
        out_name='All Data.xlsx',
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
