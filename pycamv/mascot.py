"""
Provides functionality for interacting with MASCOT data.
"""

import os
import re
import xml.etree.ElementTree as ET

from scipy.misc import comb

from pyproteome import paths
from . import ms_labels


MASCOT_NS = {
    "mascot":
        "http://www.matrixscience.com/xmlns/schema/mascot_search_results_2",
}


class PeptideQuery:
    """
    Attributes
    ----------
    gi : str
    protein : str
    query : int
    pep_rank : int
    pep_score : float
    pep_exp_mz : float
    pep_exp_z : int
    pep_seq : str
    pep_var_mods : list of tuple of (int, str, tuple of str)
    pep_fixed_mods : list of tuple of (int, str, tuple of str)
    scan : int
    num_comb : int
    """
    def __init__(
        self, gi, protein, query,
        pep_rank, pep_score, pep_exp_mz, pep_exp_z, pep_seq,
        pep_var_mods, pep_fixed_mods, scan,
    ):
        self.gi = gi
        self.protein = protein
        self.query = query
        self.pep_rank = pep_rank
        self.pep_score = pep_score
        self.pep_exp_mz = pep_exp_mz
        self.pep_exp_z = pep_exp_z
        self.pep_seq = pep_seq
        self.pep_var_mods = pep_var_mods
        self.pep_fixed_mods = pep_fixed_mods
        self.scan = scan
        self.num_comb = self._calc_num_comb()

    def _unique_tuple(self):
        return (
            self.gi,
            self.query,
            self.pep_seq,
            self.scan,
        )

    def __hash__(self):
        return hash(self._unique_tuple())

    def __eq__(self, other):
        if not isinstance(other, PeptideQuery):
            raise TypeError(other)
        return self._unique_tuple() == other._unique_tuple()

    @property
    def pep_mods(self):
        return self.pep_var_mods + self.pep_fixed_mods

    @property
    def get_label_mods(self):
        return [
            mod
            for _, mod, letters in self.pep_mods
            if mod in ms_labels.LABEL_NAMES and "N-term" in letters
        ]

    def _calc_num_comb(self):
        num_comb = 1

        for count, mod, letters in self.pep_var_mods:
            if mod == "Phospho" and letters == ["S", "T"]:
                letters = ["S", "T", "Y"]

            potential_mod_sites = sum(self.pep_seq.count(i) for i in letters)

            # Subtract sites that will be taken up by another modification
            # (i.e. Oxidation and Dioxidation of M)
            for o_count, o_mod, o_letters in self.pep_var_mods:
                if o_letters == letters and o_mod != mod:
                    potential_mod_sites -= o_count

            num_comb *= comb(potential_mod_sites, count)

        return num_comb


def _count_instances(pep_seq, letters):
    return sum(
        (["N-term"] + list(pep_seq) + ["C-term"]).count(letter)
        for letter in letters
    )


def _parse_letters(letters):
    """
    Turns a string of residue letters (i.e. "STY") into a list.

    Includes special provisions for N- and C-term modifications.

    Parameters
    ----------
    letters : str

    Returns
    -------
    tuple of str
    """
    if letters in ["N-term", "C-term"]:
        return (letters,)

    return tuple(letters)


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
            query = int(peptide.get("query"))
            rank = int(peptide.get("rank"))
            pep_score = float(peptide.find("mascot:pep_score", MASCOT_NS).text)
            exp_mz = float(peptide.find("mascot:pep_exp_mz", MASCOT_NS).text)
            exp_z = int(peptide.find("mascot:pep_exp_z", MASCOT_NS).text)
            pep_seq = peptide.find("mascot:pep_seq", MASCOT_NS).text

            var_mods = peptide.find("mascot:pep_var_mod", MASCOT_NS).text

            if fixed_mods:
                pep_fixed_mods = [
                    re.match(
                        "(.+) \((.+)\)",
                        mod.strip()
                    ).group(1, 2)
                    for mod in fixed_mods
                ]
                pep_fixed_mods = [
                    (name, _parse_letters(letters))
                    for name, letters in pep_fixed_mods
                ]
                pep_fixed_mods = [
                    (_count_instances(pep_seq, letters), name, letters)
                    for name, letters in pep_fixed_mods
                ]
                pep_fixed_mods = [
                    (count, name, letters)
                    for count, name, letters in pep_fixed_mods
                    if count > 0
                ]
            else:
                pep_fixed_mods = []

            if var_mods:
                # i.e. "2 Phospho (STY)""
                var_mods = [
                    re.match(
                        "((\d+) )?(.+) \((.+)\)",
                        mod.strip()
                    ).group(2, 3, 4)
                    for mod in var_mods.split(";")
                ]
                var_mods = [
                    (int(count) if count else 1, name, _parse_letters(letters))
                    for count, name, letters in var_mods
                ]
            else:
                var_mods = []

            scan = int(
                re.search(
                    "(scans:|Cmpd_)(\d+)",
                    peptide.find("mascot:pep_scan_title", MASCOT_NS).text
                ).group(2)
            )

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
                    pep_fixed_mods,
                    scan,
                )
            )
            index += 1

    out = sorted(out, key=lambda x: x.scan)

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
    out : list of :class:`PeptideQuery<pycamv.mascot.PeptideQuery>`
    """
    xml_path = os.path.join(paths.MASCOT_XML_DIR, xml_name)
    tree = ET.parse(xml_path)
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
