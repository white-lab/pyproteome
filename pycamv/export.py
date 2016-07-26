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


def _extract_pep_seq(sequence):
    return "".join(
        letter
        for letter, _ in sequence
        if letter not in ["N-term", "C-term"]
    )

def _extract_mods(sequence):
    return tuple(
        tuple(mods)
        for _, mods in sequence
    )


def _count_mods(sequence):
    return sum(
        bool(mods)
        for letter, mods in sequence
        if letter not in ["N-term", "C-term"]
    )


def _get_mods_description(mods):
    return "+{}".format(_count_mods(mods))


def _mod_positions(sequence):
    return [
        index
        for index, (letter, _) in enumerate(sequence)
        if letter not in ["N-term", "C-term"]
    ]


def _pep_mod_name(sequence):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in sequence
        if letter not in ["N-term", "C-term"]
    )


def _get_labels_mz(query):
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
    ### Mappings between proteins, peptides, modifications, queries, and peaks
    # Mapping for protein -> peptides
    prot_dict = DefaultOrderedDict(list)

    for (query, _), _ in peak_hits:
        prot_dict[query.protein].append(query.pep_seq)

    # Mapping for peptides -> modifications
    pep_dict = DefaultOrderedDict(list)

    for (query, sequence), _ in peak_hits.items():
        pep_dict[query.pep_seq].append(_extract_mods(sequence))

    # Mapping for modifications -> queries / sequences
    mods_dict = {
        (seq, _extract_mods(sequence): (query, sequence)
        for (query, seq), _ in peak_hits.items()
    }

    # Mapping for queries -> sequence + modifications
    query_dict = DefaultOrderedDict(list)

    for (query, seq), hits in peak_hits.items():
        query_dict[query, seq].append(_extract_mods(sequence))

    # Mapping for queries -> peak hits
    scan_data = {
        query: hits
        for (query, _), hits in peak_hits.items()
    }

    ### Pre-calculate IDs for later reference
    # Protein IDs
    prot_index = {}
    for index, prot_name in enumerate(prot_dict.keys()):
        prot_index[prot_name] = index

    # Peptide IDs
    pep_index = {}
    for index, pep_seq in enumerate(pep_dict.keys()):
        pep_index[pep_seq] = index

    # Peptide Data IDs
    # Mod State IDs
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

    # Modification IDs
    mod_index = {}

    for pep_seq, mods in pep_dict.items():
        for index, mod in enumerate(mods):
            mod_index[pep_seq, mod] = index

    scan_index = {}

    for pep_seq, mods in pep_dict.items():
        for mod in mods:
            for index, (query, _) in enumerate(mods_dict[mod]):
                scan_index[query] = index

    # Peak Match IDs
    match_index = {}

    for (_, seq), hits in peak_hits.items():
        index = 0

        for peak_hit in hits:
            for name in peak_hit.match_list.keys():
                match_index[seq, name] = index
                index += 1

    ### Individual data parsing functions
    def _gen_match_data(seq, peaks):
        """
        Generate a list of all potential peak matches
        """
        for peak_hit in peaks:
            for name, (mz, _) in peak_hit.match_list.items():
                name_split = name.split("_")
                ion_type, ion_pos = None, None

                if name_split[0] in "abc":
                    ion_type, ion_pos = "b", int(name_split[1])
                elif name_split[0] in "xyz":
                    ion_type, ion_pos = "y", int(name_split[1])

                yield OrderedDict([
                    ("id", match_index[seq, name]),
                    ("mz", mz),
                    ("name", name),
                    ("ionType", ion_type),
                    ("ionPosition", ion_pos),
                ])

    def _get_match_data(seq, peaks):
        return list(_gen_match_data(seq, peaks))

    def _get_mod_data(pep_seq, mod):
        return OrderedDict([
            ("id", mod_index[pep_seq, mod]),
            ("position", _mod_positions(mods_dict[mod][1])),
            ("name", _pep_mod_name(mods_dict[mod][1])),
            (
                "matchData",
                _get_match_data(pep_seq, scan_data[mods_dict[mod][0]]),
            ),
        ])

    def _get_mod_states(pep_seq, mods):
        return [
            OrderedDict([
                ("id", mod_index[pep_seq, mod]),
                ("modDesc", _get_mods_description(mods_dict[mod][1])),
                ("mods", _get_mod_data(mod)),
            ])
            for mod in mods
        ]

    def _get_peptide_data():
        """
        Return all information mapping (modified) peptides to their sequence,
        descriptions, ion fragmentation patterns, etc.
        """
        return [
            OrderedDict([
                ("id", pep_index[pep_seq])
                ("peptideSequence", pep_seq),
                ("modificationStates", _get_mod_states(pep_seq, mods)),
            ])
            for pep_seq, mods in pep_dict.items()
        ]

    #
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
                    "matchInfo",
                    [
                        OrderedDict([
                            ("modsId", mod_index[seq, mod]),
                            (
                                "matchId",
                                match_index.get(
                                    (seq, peak_hit.name),
                                    None,
                                )
                        ])
                        for mod in query_dict[query, seq]
                    ],
                ),
            ])
            for peak_hit in peak_hits[query, seq]
        ]

    def _get_scans(peptide):
        """
        Return information on individual scans, including peaks, precursor ions,
        and peptide modification assignments.
        """
        return [
            OrderedDict([
                ("scanNumber", query.scan),
                ("scanId", scan_index[query]),
                ("chargeState", query.pep_exp_z),
                (
                    "precursorScanData",
                    _peaks_to_dict(precursor_windows[query]),
                ),
                ("precursorMz", query.pep_exp_mz),
                ("quantScanData", _peaks_to_dict(label_windows[query])),
                ("quantMz", _get_labels_mz(query)),
                ("choiceData", _get_default_choice_data(peptide)),
                ("scanData", _get_scan_assignments(query)),
            ])
            for query in mods_dict[_extract_mods(sequence)]
        ]

    def _get_peptide_scan_data(peptides):
        """
        Map peptides to their data IDs, scans, and candidate modification
        states.
        """
        return [
            OrderedDict([
                ("peptideId", pep_index[peptide.pep_seq]),
                ("peptideDataId", pep_data_index[peptide]),
                ("modificationStateId", mod_state_index[peptide]),
                ("scans", _get_scans(peptide)),
            ])
            for peptide in peptides
        ]

    #
    def _get_scan_data():
        """
        Return all information mapping proteins / peptides to their scans and
        a list of candidate modification patterns.
        """
        return [
            OrderedDict([
                ("proteinName", prot_name),
                ("proteinId", prot_index[prot_name]),
                ("peptides", _get_peptide_scan_data(peptides)),
            ])
            for prot_name, peptides in prot_dict.items()
        ]

    data = OrderedDict([
        ("peptideData", _get_peptide_data()),
        ("scanData", _get_scan_data()),
    ])

    with open(out_path, "w") as f:
        json.dump(data, f)

    return data
