"""
Provides functionality for exporting processed data to .camv files (JSON
format).
"""

from collections import OrderedDict
import json
import re

from .utils import DefaultOrderedDict

from . import ms_labels
from pyproteome.loading import RE_PROTEIN


RE_B_Y_IONS = re.compile("([abcxyz]_\{[0-9]+\})\^\{\+\}")
SUPERSCRIPT_UNICODE_START = ord("\u2070")
SUBSCRIPT_UNICODE_START = ord('\u2080')
SCRIPT_MAPPING = {
    str(i): i
    for i in range(10)
}
SCRIPT_MAPPING["+"] = 10
SCRIPT_MAPPING["-"] = 11
SCRIPT_MAPPING["("] = 12
SCRIPT_MAPPING[")"] = 13


def _rewrite_ion_name(name):
    m = RE_B_Y_IONS.match(name)

    if m:
        name = m.group(1)

    ret = ""
    sup, sub = False, False
    paren = 0

    for char in name:
        if char in "^_{}":
            if char == "^":
                sup, sub = True, False
            elif char == "_":
                sup, sub = False, True
            elif char == "}":
                sup, sub = False, False
                paren -= 1
            elif char == "{":
                paren += 1
            continue

        if sup:
            if char == "1":
                ret += u"\u00B9"
            elif char == "2":
                ret += u"\u00B2"
            elif char == "3":
                ret += u"\u00B3"
            else:
                ret += chr(SUPERSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        elif sub:
            ret += chr(SUBSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        else:
            ret += char

        if sup or sub:
            if not paren:
                sup, sub = False, False

    return ret


def _peaks_to_dict(peaks):
    return [
        OrderedDict([
            ("mz", mz),
            ("into", i),
        ])
        for mz, i in peaks
    ]


def _join_seq_mods(seq, mods):
    return tuple(zip(["N-term"] + list(seq) + ["C-term"], mods))


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


def _get_mods_description(pep_seq, mod_state):
    return ("+ " if mod_state else "") + " - ".join(
        "{} {}{}".format(count, mod[0].lower(), "".join(letters))
        for count, mod, letters in mod_state
    )


def _mod_positions(mods):
    return [
        index
        for index, mod in enumerate(mods)
        if mod
    ]


def _pep_mod_name(pep_seq, mods):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in zip(pep_seq, mods[1:-1])
    )


def _get_labels_mz(query):
    return [
        mz
        for mod in set(query.get_label_mods)
        for mz in ms_labels.LABEL_MASSES.get(mod, [])
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

    Example
    -------
    >>> from pycamv import export, validate
    >> name = "2016-06-27-CKX13-pY-NTA-elute-pre67-col77"
    >>> options, peak_hits, precursor_windows, label_windows = \
    ...     validate.validate_spectra(name)
    >>> json_out = export.export_to_camv(
    ...     "{}.camv".format(name),
    ...     peak_hits, precursor_windows, label_windows,
    ... )
    """
    ###
    # Mappings between proteins, peptides, modifications, queries, and peaks
    ###
    # Mapping for protein -> peptides
    prot_dict = DefaultOrderedDict(set)

    for query, _ in peak_hits.keys():
        prot_dict[RE_PROTEIN.match(query.protein).group(1)].add(query.pep_seq)

    pep_dict = DefaultOrderedDict(set)

    for query, _ in peak_hits.keys():
        pep_dict[query.pep_seq].add(tuple(query.pep_var_mods))

    mod_states_dict = DefaultOrderedDict(set)

    for query, seq in peak_hits.keys():
        mod_states_dict[query.pep_seq, tuple(query.pep_var_mods)].add(
            _extract_mods(seq)
        )

    mod_query_dict = DefaultOrderedDict(set)

    for query, seq in peak_hits.keys():
        mod_query_dict[query.pep_seq, tuple(query.pep_var_mods)].add(query)

    # Mapping for modifications -> queries
    mods_dict = DefaultOrderedDict(list)

    for query, seq in peak_hits.keys():
        mods_dict[query.pep_seq, _extract_mods(seq)].append(query)

    # Mapping for queries -> sequence + modifications
    query_dict = DefaultOrderedDict(list)

    for query, seq in peak_hits.keys():
        query_dict[query].append((query.pep_seq, _extract_mods(seq)))

    ###
    # Pre-calculate IDs for later reference
    ###
    # Protein IDs
    prot_index = {
        prot_name: index
        for index, prot_name in enumerate(sorted(prot_dict.keys()))
    }

    # Peptide IDs
    pep_index = {
        (pep_seq, mod_states): index
        for peptides in prot_dict.values()
        for index, (pep_seq, mod_states) in enumerate(
            (pep_seq, mod_states)
            for pep_seq in peptides
            for mod_states in pep_dict[pep_seq]
        )
    }

    # Peptide Data IDs
    pep_data_index = {
        pep_seq: index
        for index, pep_seq in enumerate(
            sorted(
                pep_seq
                for peptides in prot_dict.values()
                for pep_seq in peptides
            )
        )
    }

    # Mod State IDs
    mod_state_index = {
        (pep_seq, mod_state): index
        for pep_seq, mod_states in pep_dict.items()
        for index, mod_state in enumerate(mod_states)
    }

    # Modification IDs
    mod_index = {
        (pep_seq, mod): index
        for (pep_seq, _), mods in mod_states_dict.items()
        for index, mod in enumerate(mods)
    }

    # Scan IDs
    scan_index = {
        query: index
        for pep_seq, mod_states in pep_dict.items()
        for mod_state in mod_states
        for index, query in enumerate(
            sorted(
                mod_query_dict[pep_seq, mod_state],
                key=lambda x: x.scan,
            )
        )
    }

    # Peak Match IDs
    match_index = {
        (pep_seq, mods, name): index
        for (pep_seq, mods), queries in mods_dict.items()
        for index, name in enumerate(
            name
            for name in sorted(
                set(
                    name
                    for query in queries
                    for peak_hit in peak_hits[
                        query, _join_seq_mods(pep_seq, mods)
                    ]
                    if peak_hit.match_list
                    for name in peak_hit.match_list.keys()
                )
            )
        )
    }

    ###
    # Individual data parsing functions
    ###
    def _gen_match_data(pep_seq, mods):
        """
        Generate a list of all potential peak matches
        """
        visited = set()
        queries = mods_dict[pep_seq, mods]

        for query in queries:
            for peak_hit in peak_hits[query, _join_seq_mods(pep_seq, mods)]:
                if not peak_hit.match_list:
                    continue

                for name, (mz, _) in peak_hit.match_list.items():
                    if name in visited:
                        continue

                    name_split = re.split("[_\^]", name)
                    name_split = [i.strip("{}") for i in name_split]

                    ion_type, ion_pos = None, None

                    if len(name_split) > 1:
                        if name_split[0] in "abc":
                            ion_type, ion_pos = "b", int(name_split[1])
                        elif name_split[0] in "xyz":
                            ion_type, ion_pos = "y", int(name_split[1])

                    yield OrderedDict([
                        ("id", match_index[pep_seq, mods, name]),
                        ("mz", mz),
                        ("name", _rewrite_ion_name(name)),
                        ("ionType", ion_type),
                        ("ionPosition", ion_pos),
                    ])

                    visited.add(name)

    def _get_match_data(pep_seq, mods):
        return sorted(
            _gen_match_data(pep_seq, mods),
            key=lambda x: x["id"],
        )

    def _get_mod_data(pep_seq, mods):
        return [
            OrderedDict([
                ("id", mod_index[pep_seq, mod]),
                ("position", _mod_positions(mod)),
                ("name", _pep_mod_name(pep_seq, mod)),
                ("matchData", _get_match_data(pep_seq, mod)),
            ])
            for mod in sorted(
                mods,
                key=lambda x: mod_index[pep_seq, x],
            )
        ]

    def _get_mod_states(pep_seq, mod_states):
        return [
            OrderedDict([
                ("id", mod_state_index[pep_seq, mod_state]),
                ("modDesc", _get_mods_description(pep_seq, mod_state)),
                (
                    "mods",
                    _get_mod_data(
                        pep_seq,
                        mod_states_dict[pep_seq, mod_state],
                    ),
                ),
            ])
            for mod_state in sorted(
                mod_states,
                key=lambda x: mod_state_index[pep_seq, x],
            )
        ]

    def _get_peptide_data():
        """
        Return all information mapping (modified) peptides to their sequence,
        descriptions, ion fragmentation patterns, etc.
        """
        return [
            OrderedDict([
                ("id", pep_data_index[pep_seq]),
                ("peptideSequence", pep_seq),
                ("modificationStates", _get_mod_states(pep_seq, mod_states)),
            ])
            for pep_seq, mod_states in sorted(
                pep_dict.items(),
                key=lambda x: pep_data_index[x[0]],
            )
        ]

    def _get_default_choice_data(pep_seq, mod_state):
        return [
            OrderedDict([
                ("modsId", mod_index[pep_seq, mod]),
                ("state", None),  # null
            ])
            for mod in sorted(
                mod_states_dict[pep_seq, mod_state],
                key=lambda x: mod_index[pep_seq, x],
            )
        ]

    def _get_matches(peak_index, query):
        return [
            OrderedDict([
                ("modsId", mod_index[seq, mods]),
                (
                    "matchId",
                    match_index.get(
                        (
                            seq,
                            mods,
                            peak_hits[
                                query,
                                _join_seq_mods(seq, mods),
                            ][peak_index].name,
                        ),
                        None,
                    ),
                ),
            ])
            for seq, mods in sorted(
                query_dict[query],
                key=lambda x: mod_index[x],
            )
        ]

    def _get_scan_assignments(query, seq):
        mod = query_dict[query][0][1]

        return [
            OrderedDict([
                ("mz", peak_hit.mz),
                ("into", peak_hit.intensity),
                ("matchInfo", _get_matches(peak_index, query)),
            ])
            for peak_index, peak_hit in sorted(
                enumerate(peak_hits[query, _join_seq_mods(seq, mod)]),
                key=lambda x: x[1].mz,
            )
        ]

    def _get_scans(pep_seq, mod_state):
        """
        Return information on individual scans, including peaks, precursor
        ions, and peptide modification assignments.
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
                (
                    "choiceData",
                    _get_default_choice_data(pep_seq, mod_state),
                ),
                (
                    "scanData",
                    _get_scan_assignments(query, pep_seq),
                ),
            ])
            for query in sorted(
                mod_query_dict[pep_seq, mod_state],
                key=lambda x: scan_index[x],
            )
        ]

    def _get_peptide_scan_data(peptides):
        """
        Map peptides to their data IDs, scans, and candidate modification
        states.
        """
        return [
            OrderedDict([
                ("peptideId", pep_index[pep_seq, mod_state]),
                ("peptideDataId", pep_data_index[pep_seq]),
                ("modificationStateId", mod_state_index[pep_seq, mod_state]),
                ("scans", _get_scans(pep_seq, mod_state)),
            ])
            for pep_seq, mod_state in sorted(
                (
                    (pep_seq, mod_state)
                    for pep_seq in peptides
                    for mod_state in pep_dict[pep_seq]
                ),
                key=lambda x: pep_index[x],
            )
        ]

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
            for prot_name, peptides in sorted(
                prot_dict.items(),
                key=lambda x: prot_index[x[0]],
            )
        ]

    data = OrderedDict([
        ("peptideData", _get_peptide_data()),
        ("scanData", _get_scan_data()),
    ])

    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)

    return data
