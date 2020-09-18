'''
This module provides functionality for interfacing with protein data.
'''

import logging

import pyproteome as pyp


LOGGER = logging.getLogger('pyproteome.protein')


class Proteins:
    '''
    Wraps a list of proteins.

    Attributes
    ----------
    proteins : tuple of :class:`.Protein`
        List of proteins to which a peptide sequence is mapped.
    '''

    def __init__(self, proteins=None):
        if proteins is None:
            proteins = ()

        self.proteins = tuple(sorted(proteins))

    def __iter__(self):
        return iter(self.proteins)

    def __len__(self):
        return len(self.proteins)

    def __hash__(self):
        return hash(
            self.proteins,
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

    def __lt__(self, other):
        return self.proteins < other.proteins

    def __str__(self):
        return ' / '.join(
            str(i)
            for i in self.proteins
        )

    @property
    def accessions(self):
        '''
        List of UniPort accessions for a group of proteins.

        Returns
        -------
        tuple of str
        '''
        return tuple(i.accession for i in self.proteins)

    @property
    def descriptions(self):
        '''
        List of protein descriptions for a group of proteins.

        Returns
        -------
        tuple of str
        '''
        return tuple(i.description for i in self.proteins)

    @property
    def genes(self):
        '''
        List of UniPort gene names for a group of proteins.

        Returns
        -------
        tuple of str
        '''
        return tuple(i.gene for i in self.proteins)


class Protein:
    '''
    Contains information about a single protein.

    Attributes
    ----------
    accession : str
        The UniProt accession (i.e. 'P40763').
    gene : str
        The UniProt gene name (i.e. 'STAT3').
    description : str
        A brief description of the protein (i.e. 'Signal transducer and
        activator of transcription 3').
    full_sequence : str
        The full sequence of the protein.
    '''

    def __init__(
        self,
        accession=None,
        gene=None, 
        description=None,
        full_sequence=None,
    ):
        self.accession = accession
        self.gene = gene
        self.description = description
        self.full_sequence = full_sequence

        if any(i is None for i in [gene, description, full_sequence]):
            up_data = pyp.pypuniprot.get_uniprot_data(accession)

            if 'gene' in up_data:
                self.gene = up_data['gene']
            elif 'id' in up_data:
                self.gene = up_data['id']
            else:
                LOGGER.warning(
                    'Unable to find {} in uniprot db'.format(accession)
                )
                self.gene = accession

            self.description = up_data.get('descriptions', [''])[0]
            self.full_sequence = up_data.get('sequence', None)

    def __hash__(self):
        return hash(self.accession)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            raise TypeError()

        return self.accession == other.accession

    def __lt__(self, other):
        return self.gene < other.gene

    def __str__(self):
        return '{} ({})'.format(
            self.description,
            self.gene,
        )
