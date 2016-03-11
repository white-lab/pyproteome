"""
Provides functionality for interacting with MASCOT data.
"""

import os
from collections import namedtuple
import re
import xml.etree.ElementTree as ET

from pyproteome import paths

MASCOT_NS = {
    "mascot":
        "http://www.matrixscience.com/xmlns/schema/mascot_search_results_2",
}

PeptideQuery = namedtuple(
    "PeptideQuery",
    [
        "gi",
        "protein",
        "query",
        "pep_rank",
        "pep_score",
        "pep_exp_mz",
        "pep_exp_z",
        "pep_seq",
        "pep_var_mods",
        "scan",
    ],
)


def _parse_mascot_2_4_1(root):
    fixed_mods = [
        i.text
        for i in root.findall(
            "mascot:fixed_mods/mascot:modification/mascot:name",
            MASCOT_NS,
        )
    ]
    variable_mods = [
        i.text
        for i in root.findall(
            "mascot:variable_mods/mascot:modification/mascot:name",
            MASCOT_NS,
        )
    ]

    scan_used = {}
    index = 0
    out = []

    for hit in root.findall("mascot:hits/mascot:hit", MASCOT_NS):
        accession = hit.find("mascot:protein", MASCOT_NS).get("accession")
        prot_desc = hit.find("mascot:protein/mascot:prot_desc", MASCOT_NS).text

        for peptide in hit.findall("mascot:protein/mascot:peptide", MASCOT_NS):
            query = peptide.get("query")
            rank = peptide.get("rank")
            pep_score = peptide.find("mascot:pep_score", MASCOT_NS).text
            exp_mz = peptide.find("mascot:pep_exp_mz", MASCOT_NS).text
            exp_z = peptide.find("mascot:pep_exp_z", MASCOT_NS).text
            pep_seq = peptide.find("mascot:pep_seq", MASCOT_NS).text

            var_mods = peptide.find("mascot:pep_var_mod", MASCOT_NS).text

            if var_mods:
                var_mods = [
                    re.match(
                        "((\d+) )?(.*)",
                        mod.strip()
                    ).group(2, 3)
                    for mod in var_mods.split(";")
                ]
                var_mods = [
                    (count if count else 1, name)
                    for count, name in var_mods
                ]
            else:
                var_mods = []

            scan = re.search(
                "(scans:|Cmpd_)(\d+)",
                peptide.find("mascot:pep_scan_title", MASCOT_NS).text
            ).group(2)

            # scan_used
            if scan in scan_used:
                index_match, rank_match = scan_used[scan]

                if rank >= rank_match:
                    continue
                else:
                    del out[index_match]

            scan_used[scan] = (index, rank)
            out.append(
                PeptideQuery(
                    accession,
                    prot_desc,
                    query,
                    rank,
                    pep_score,
                    exp_mz,
                    exp_z,
                    pep_seq,
                    var_mods,
                    scan,
                )
            )
            index += 1

    return fixed_mods, variable_mods, out


def read_mascot_xml(xml_name):
    """
    Parse a MASCOT XML file.

    Parameters
    ----------
    xml_name : str
        Name of XML file in Mascot XMLs directory.

    Returns
    -------
    fixed_mods : list of str
    var_mods : list of str
    out : list of pyproteome.mascot.PeptideQuery
    """
    tree = ET.parse(os.path.join(paths.MASCOT_XML_DIR, xml_name))
    root = tree.getroot()

    ms_version = root.find("mascot:header/mascot:MascotVer", MASCOT_NS).text
    ms_version = tuple(int(i) for i in ms_version.split("."))

    if ms_version >= (2, 4, 1):
        return _parse_mascot_2_4_1(root)
#     elif ms_version == "2.4.0":
#         return _parse_mascot_2_4_0(root)
#     elif ms_version == "2.3.02":
#         return _parse_mascot_2_3_02(root)
#     elif ms_version == "2.1.03":
#         return _parse_mascot_2_1_03(root)
    else:
        raise Exception("Unsupported Mascot Version: {}".format(ms_version))
