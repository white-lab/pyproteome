"""
Provides functionality for exporting processed data to .camv files (JSON
format).
"""

from collections import defaultdict, OrderedDict
import json


def _peaks_to_dict(peaks):
    return [
        OrderedDict([
            ("mz", mz),
            ("into", i),
        ])
        for mz, i in peaks
    ]


def _extract_mods(sequence):
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
    prot_dict = defaultdict(list)
    for (query, _), _ in peak_hits:
        prot_dict[query.protein].append(query.pep_seq)

    # Mapping for peptides -> modifications
    pep_dict = defaultdict(list)
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

    def _get_match_data(peaks):
        match_index = 0

        # XXX: ion_position start?
        for ion_position, peak_hit in enumerate(peaks, start=1):
            if not peak_hit.name:
                continue

            yield OrderedDict([
                ("id", match_index),
                ("mz", peak_hit.mz),
                ("name", peak_hit.name),
                ("ionType", "b")  # XXX: FIXME!
                ("ionPosition", ion_position),
            ])

            match_index += 1

    def _mod_states(mod_index, mod):
        return OrderedDict([
            ("id", mod_index),
            ("modDesc", "+{}".format(_count_mods(mods_dict[mod][1]))),
            (
                "mods",
                OrderedDict([
                    ("id", None),
                    ("position", _mod_positions(mods_dict[mod][1])),
                    ("name", _mod_name(mods_dict[mod][1])),
                    (
                        "matchData",
                        list(_get_match_data(scan_data[mods_dict[mod][0]])),
                    ),
                ]),
            ),
        ])

    data["peptideData"] = [
        OrderedDict([
            ("id", pep_index)
            ("peptideSequence", pep_seq),
            (
                "modificationStates",
                [
                    _mod_states(mod_index, mod)
                    for mod_index, mod in enumerate(mods)
                ],
            ),
        ])
        for pep_index, (pep_seq, mods) in enumerate(pep_dict.items())
    ]

    data["scanData"] = [
        OrderedDict([
            ("proteinName", protein_name),
            ("proteinId", prot_index),
            (
                "peptides",
                [
                    OrderedDict([
                        ("peptideId", query_index),
                        ("peptideDataId", None),
                        ("modificationStateId", None),
                        (
                            "scans",
                            [
                                OrderedDict([
                                    ("scanNumber", None),
                                    ("scanId", None),
                                    ("chargeState", None),
                                    (
                                        "precursorScanData",
                                        _peaks_to_dict(
                                            precursor_windows[query]
                                        )
                                    ),
                                    ("precursorMz", None),
                                    (
                                        "quantScanData",
                                        _peaks_to_dict(
                                            label_windows[query]
                                        )
                                    ),
                                    ("quantMz", None),
                                    (
                                        "choiceData",
                                        [
                                            OrderedDict([
                                                ("modsId", mod_index),
                                                ("state", None),  # null
                                            ])
                                            for mod_index, _ in mods_dict[query]
                                        ]
                                    ),
                                    (
                                        "scanData",
                                        [
                                            OrderedDict([
                                                ("mz", None),
                                                ("into", None),
                                                (
                                                    "matchInfo",
                                                    OrderedDict([
                                                        ("modsId", None),
                                                        ("matchId", None),
                                                    ]),
                                                ),
                                            ])
                                        ],
                                    ),
                                ]),
                            ]
                        ),
                    ])
                    for pep_index, peptide in peptides
                ]
            ),
        ])
        for prot_index, (protein_name, peptides) in enumerate(prot_dict.items())
    ]

    with open(out_path, "w") as f:
        json.dump(data, f)

    return data
