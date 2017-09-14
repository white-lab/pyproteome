import os

try:
    import Bio.Alphabet.IUPAC
    import Bio.Seq
    import Bio.motifs
except ImportError:
    Bio = None

from . import motif


def make_logo(
    data, f,
    folder_name=None, filename="Motif.svg",
    **kwargs
):
    """
    Create a sequence logo figure.

    Logos are created based on the frequencies of peptides in a data set.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    filename : str, optional
    """
    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    nmer_args = motif.get_nmer_args(kwargs)

    alpha = Bio.Alphabet.IUPAC.protein
    m = Bio.motifs.create(
        [
            Bio.Seq.Seq(seq.upper(), alphabet=alpha)
            for seq in motif.generate_n_mers(
                data.filter(**f).psms["Sequence"],
                **nmer_args
            )
        ],
        alphabet=alpha,
    )
    m.weblogo(
        filename,
        yaxis_scale=3,
        format="svg",
        # logo_title="Testing",
        stack_width="large"
    )
