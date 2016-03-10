"""
This module tracks the path to user data files. Developers can override paths
here when using a custom data hierarchy.
"""

import os

BASE_DIR = os.path.abspath("..")

BCA_ASSAY_DIR = os.path.join(BASE_DIR, "BCA Protein Assays")
CAMV_OUT_DIR = os.path.join(BASE_DIR, "CAMV Output")
CAMV_SESS_DIR = os.path.join(BASE_DIR, "CAMV Sessions")
MASCOT_XML_DIR = os.path.join(BASE_DIR, "Mascot XMLs")
MS_SEARCHED_DIR = os.path.join(BASE_DIR, "MS Searched")
MS_RAW_DIR = os.path.join(BASE_DIR, "MS RAW")
SCAN_LISTS_DIR = os.path.join(BASE_DIR, "Scan Lists")
SCRIPTS_DIR = os.path.join(BASE_DIR, "Scripts")
