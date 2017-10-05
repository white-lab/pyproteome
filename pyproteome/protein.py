"""
This module provides functionality for interfacing with protein data.
"""

import logging

from . import fetch_data


LOGGER = logging.getLogger("pyproteome.protein")


class Proteins:
    """
    Wraps a list of proteins.

    Attributes
    ----------
    proteins : list of :class:`Protein<pyproteome.protein.Protein>`
    """

    def __init__(self, proteins):
        self.proteins = proteins

    def __iter__(self):
        return iter(self.proteins)

    def __hash__(self):
        return hash(
            tuple(self.proteins),
        )

    def __eq__(self, other):
        if isinstance(other, str):
            return (
                any(i == other for i in self.genes) or
                any(i == other for i in self.accessions)
            )

        if not isinstance(other, Proteins):
            raise TypeError(type(other))

        return len(self.proteins) == len(other.proteins) and all(
            i == j
            for i, j in zip(self.proteins, other.proteins)
        )

    def __str__(self):
        return " / ".join(
            str(i)
            for i in self.proteins
        )

    @property
    def accessions(self):
        return [i.accession for i in self.proteins]

    @property
    def genes(self):
        return [i.gene for i in self.proteins]


class Protein:
    """
    Contains information about a single protein.

    Attributes
    ----------
    accession : str
    gene : str
    description : str
    full_sequence : str
    """

    def __init__(
        self, accession,
        gene=None, description=None, full_sequence=None
    ):
        self.accession = accession
        self.gene = gene
        self.description = description
        self.full_sequence = full_sequence

        if any(i is None for i in [gene, description, full_sequence]):
            up_data = fetch_data.get_uniprot_data(accession)

            if "gene" in up_data:
                self.gene = up_data["gene"]
            elif "id" in up_data:
                self.gene = up_data["id"]
            else:
                LOGGER.warning(
                    "Unable to find {} in uniprot db".format(accession)
                )
                self.gene = accession

            self.description = up_data.get("descriptions", [""])[0]
            self.full_sequence = up_data.get("sequence", None)

    def __hash__(self):
        return hash(self.accession)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            raise TypeError()

        return self.accession == other.accession

    def __str__(self):
        return "{} ({})".format(
            self.description,
            self.gene,
        )
