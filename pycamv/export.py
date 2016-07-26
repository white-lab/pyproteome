"""
Provides functionality for exporting processed data to .camv files (JSON
format).
"""

from collections import OrderedDict
import json

from .utils import DefaultOrderedDict

from . import ms_labels


def _peaks_to_dict(peaks):
    return [
        OrderedDict([
            ("mz", mz),
            ("into", i),
        ])
        for mz, i in peaks
    ]


def _extract_mods(sequence):
    # XXX: Fill me in!
    return ""


def _count_mods(sequence):
    return sum(
        bool(mods)
        for letter, mods in sequence
        if letter not in ["N-term", "C-term"]
    )


def _mod_positions(sequence):
    return [
        index
        for index, (letter, mods) in enumerate(sequence)
        if letter not in ["N-term", "C-term"]
    ]


def _mod_name(sequence):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in sequence
        if letter not in ["N-term", "C-term"]
    )


def _get_label_mz(query):
    return [
        mz
        for var_mod in query.pep_var_mods
        for mz in ms_labels.LABEL_MASSES.get(var_mod, [])
    ]


def export_to_camv(out_path, peak_hits, precursor_windows, label_windows):
    """
    Export data to .camv file.

    Parameters
    ----------
    out_path : str

    Returns
    -------
    data : dict
        Dictionary of data written to file.
    """
    data = OrderedDict()

    # Mapping for protein -> peptides
    prot_dict = DefaultOrderedDict(list)
    for (query, _), _ in peak_hits:
        prot_dict[query.protein].append(query.pep_seq)

    # Mapping for peptides -> modifications
    pep_dict = DefaultOrderedDict(list)
    for (query, sequence), _ in peak_hits:
        pep_dict[query.pep_seq].append(_extract_mods(sequence))

    # Mapping for modifications -> queries / sequences
    mods_dict = {
        _extract_mods(sequence): (query, sequence)
        for (query, sequence), _ in peak_hits
    }

    # Mapping for queries -> peak hits
    scan_data = {
        query: hits
        for (query, _), hits in peak_hits
    }

    # Pre-calculate IDs for later reference
    prot_index = {}
    for index, prot_name in enumerate(prot_dict.keys()):
        prot_index[prot_name] = index

    pep_index = {}
    for index, pep_seq in enumerate(pep_dict.keys()):
        pep_index[pep_seq] = index

    pep_data_index, mod_state_index = {}, {}
    pep_index, mod_index = 0, 0

    for mods, (query, seq) in queries_dict.items():
        key = (seq.pep_seq)

        if key not in pep_data_index:
            pep_data_index[key] = pep_index
            pep_index += 1

        key = (seq.pep_seq, mods)

        if key not in mod_state_index:
            mod_state_index[key] = mod_index
            mod_index += 1

    mod_index = {}

    for index, mod in enumerate(mods_dict.keys()):
        mod_index[mod] = index

    scan_index = {}

    for pep_seq, mods in pep_dict.items():
        for mod in mods:
            for index, (query, _) in enumerate(mods_dict[mod]):
                scan_index[query] = index

    match_index = {}
    ...

    # Individual data parsing functions
    def _get_match_data(peaks):
        match_index = 0

        for peak_hit in peaks:
            if not peak_hit.name:
                continue

            ion_type, ion_pos = None, None
            name_split = peak_hit.name.split("_")

            if name_split[0] in "abc":
                ion_type, ion_pos = "b", int(name_split[1])
            elif name_split[0] in "xyz":
                ion_type, ion_pos = "y", int(name_split[1])

            yield OrderedDict([
                ("id", match_index),
                ("mz", peak_hit.mz),
                ("name", peak_hit.name),
                ("ionType", ion_type),
                ("ionPosition", ion_pos),
            ])

            match_index += 1

    def _get_mods_description(mods):
        return "+{}".format(_count_mods(mods))

    def _get_mod_data(mod):
        return OrderedDict([
            ("id", None),
            ("position", _mod_positions(mods_dict[mod][1])),
            ("name", _mod_name(mods_dict[mod][1])),
            (
                "matchData",
                list(
                    _get_match_data(scan_data[mods_dict[mod][0]])
                ),
            ),
        ])

    def _get_mod_states(mods):
        return [
            OrderedDict([
                ("id", mod_index),
                ("modDesc", _get_mods_description(mods_dict[mod][1])),
                ("mods", _get_mod_data(mod)),
            ])
            for mod_index, mod in enumerate(mods)
        ]

    def _get_default_choice_data(mods_dict):
        return [
            OrderedDict([
                ("modsId", mod_index),
                ("state", None),  # null
            ])
            for mod_index, _ in mods_dict[query]
        ]

    def _get_scan_assignments(query, seq):
        return [
            OrderedDict([
                ("mz", peak_hit.mz),
                ("into", peak_hit.intensity),
                (
                    "matchInfo", [
                        OrderedDict([
                            ("modsId", None),
                            ("matchId", None),
                        ])
                        for match in peak_hit.match_list
                    ],
                ),
            ])
            for peak_hit in peak_hits[query, seq]
        ]

    def _get_scan_data(peptide):
        return [
            OrderedDict([
                ("scanNumber", query.scan_number),
                ("scanId", scan_index[query]),
                ("chargeState", query.pep_exp_z),
                (
                    "precursorScanData",
                    _peaks_to_dict(precursor_windows[query]),
                ),
                ("precursorMz", query.pep_exp_mz),
                (
                    "quantScanData",
                    _peaks_to_dict(label_windows[query]),
                ),
                ("quantMz", _get_label_mz(query)),
                ("choiceData", _get_default_choice_data(mods_dict)),
                ("scanData", _get_scan_assignments(query)),
            ])
            for query in mods_dict[_extract_mods(sequence)]
        ]

    def _get_peptide_scan_data(peptides):
        return [
            OrderedDict([
                ("peptideId", pep_index),
                ("peptideDataId", pep_data_index[peptide]),
                ("modificationStateId", mod_state_index[peptide]),
                ("scans", _get_scan_data(peptide)),
            ])
            for pep_index, peptide in peptides
        ]

    data["peptideData"] = [
        OrderedDict([
            ("id", pep_index)
            ("peptideSequence", pep_seq),
            ("modificationStates", _get_mod_states(mods)),
        ])
        for pep_index, (pep_seq, mods) in enumerate(pep_dict.items())
    ]

    data["scanData"] = [
        OrderedDict([
            ("proteinName", prot_name),
            ("proteinId", prot_index[prot_name]),
            ("peptides", _get_peptide_scan_data(peptides)),
        ])
        for prot_name, peptides in prot_dict.items()
    ]

    with open(out_path, "w") as f:
        json.dump(data, f)

    return data
