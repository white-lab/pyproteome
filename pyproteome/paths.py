"""
This module tracks the path to user data files. Developers can override paths
here when using a custom data hierarchy.
"""

import os

BASE_DIR = os.path.abspath("..")

MS_SEARCHED_DIR = os.path.join(BASE_DIR, "MS Searched")
MS_RAW_DIR = os.path.join(BASE_DIR, "MS RAW")
MASCOT_XML_DIR = os.path.join(BASE_DIR, "Mascot XMLs")
CAMV_OUT_DIR = os.path.join(BASE_DIR, "CAMV Output")
SCAN_LISTS_DIR = os.path.join(BASE_DIR, "Scan Lists")
CAMV_SESS_DIR = os.path.join(BASE_DIR, "CAMV Sessions")
