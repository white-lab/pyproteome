from matplotlib.patches import Circle, Ellipse, Rectangle, FancyBboxPatch
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter


def _get_protein_seq(slc, gene):
    for _, row in slc.iterrows():
        for protein in row['Proteins'].proteins:
            if protein.gene == gene:
                return protein.full_sequence
            

def _get_lc(
    fc, 
    p, 
    edge=False, 
    p_cutoff=1e-2,
    upper_fc_cutoff=1.05,
    lower_fc_cutoff=0.95,
):
    if p < p_cutoff:
        if fc > upper_fc_cutoff:
            return '#ff0000' if edge else '#ff8888'
        elif fc < lower_fc_cutoff:
            return '#0000ff' if edge else '#8888ff'
    return '#000000' if edge else '#888888'


def _draw_line(
    ax, 
    row, 
    col, 
    col2, 
    larrow=True, 
    rarrow=True, 
    line_color='#000000',
    row_offset=-.5,
    col_offset=-.5,
    z=0,
):
    return ax.add_artist(
        FancyBboxPatch(
            (col + col_offset, row - .5),
            col2 - col, 
            .9,
            ec=None,
            fc=line_color,
            zorder=-1,
            boxstyle="round,pad=.05",
        )
    )


def draw_protein_seq(
    ds,
    genes,
    max_col=50,
    p_cutoff=1e-2,
    upper_fc_cutoff=1.05,
    lower_fc_cutoff=0.95,
    missed_cleavage=1,
):
    '''
    Generate a figure showing all peptides in a data set mapping to the
    full sequence of their respective proteins.
    
    Peptide differential regulation is indicated by bars for full peptide
    sequences and circles indicating phosphorylated residues.
    
    Bars and circles are colored red for upregulation and blue for downregulation.
    Dark grey bars indicate an unmodified peptide with no change. Light grey bars
    indicate that only the phosphorylated version of that peptide was identified.

    Parameters
    ----------
    ds : :class:`pyproteome.data_sets.data_set.DataSet`
    genes : list of str
    max_col : int, optional
    p_cutoff : float, optional
    upper_fc_cutoff : float, optional
    lower_fc_cutoff : float, optional
    missed_cleavage : int, optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`

    Examples
    --------
    >>> figs = analysis.protein.draw_protein_seq(
    ...     ds, ['Mapt']
    ... )    
    '''
    ds = ds.filter(protein=genes)
    ds = ds.filter(missed_cleavage=missed_cleavage)

    max_col = 50
    print(genes)

    # mod_types = sorted(set([
    #     mod.mod_type 
    #     for seq in ds['Sequence']
    #     for mod in seq.modifications 
    #     if 'TMT' not in mod.mod_type
    # ]))
    # colors = sns.color_palette(
    #     "hls", len(mod_types),
    # ).as_hex()
    mod_types = ['Phospho'] 
    colors = ['#db5f57']
    
    figs = []

    for gene in genes:
        slc = ds[ds['Proteins'] == gene]
        try:
            slc['Sort'] = slc.apply(
                lambda x: 
                (
                    len(x['Modifications'].get_mods(mod_types)) > 0,
                    len(x['Sequence']),
                    x['p-value'],
                ), 
                axis=1,
            )
            slc = slc.sort_values('Sort', ascending=False)
        except:
            continue
#         display(slc[['Proteins', 'Sequence', 'Fold Change', 'p-value']])

        prot_seq = _get_protein_seq(slc, gene)
        row_max = int(np.ceil(len(prot_seq) / max_col))

        fig_x, fig_y = row_max / 4, max_col / 5
        f, ax = plt.subplots(figsize=(fig_y, fig_x), dpi=300)
        figs.append(f)
        
        ax.set_title(gene)
        ax.set_xlim(left=-1, right=max_col)
        ax.set_ylim(bottom=row_max - .5, top=-1)
        ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
        ax.set_yticks(list(range(row_max)))
        ax.set_yticklabels([i * max_col + 1 for i in ax.get_yticks()])

        row, col = 0, 0
        for ind, letter in enumerate(prot_seq):
            ax.text(
                x=col,
                y=row,
                s=letter,
                ha='center',
                va='center',
                family='monospace',
                zorder=20,
            )
            col += 1
            if col >= max_col:
                col = 0
                row += 1

        count = Counter()

        for _, pep in slc.iterrows():
            seq = pep['Sequence']
            match = [i for i in seq.protein_matches if i.protein.gene == gene][0]

            fc_lc = _get_lc(
                pep['Fold Change'], 
                pep['p-value'], 
                p_cutoff=p_cutoff,
                upper_fc_cutoff=upper_fc_cutoff,
                lower_fc_cutoff=lower_fc_cutoff,
            )
            ec_lc = _get_lc(
                pep['Fold Change'], 
                pep['p-value'], 
                edge=True, 
                p_cutoff=p_cutoff,
                upper_fc_cutoff=upper_fc_cutoff,
                lower_fc_cutoff=lower_fc_cutoff,
            )

            mods = seq.modifications.get_mods(mod_types)
            # print(seq.pep_seq, match.rel_pos, match.exact, str(seq))

            row = match.rel_pos // max_col
            col = match.rel_pos % max_col
            col2 = (col + len(seq.pep_seq))

            count = Counter(seq.pep_seq)

            ro = -.4 - .1 * count[row, col]
            ro = -.4 - .1 * (count['K'] + count['R'])

            # Circle phosphosites
            for mod in mods:
                letter = mod.letter
                ind = mod.rel_pos + match.rel_pos
                mod_row = ind // max_col
                mod_col = ind % max_col

                ax.text(
                    x=mod_col,
                    y=mod_row,
                    s=letter,
                    ha='center',
                    va='center',
                    color='k',
                    family='monospace',
                    zorder=20,
                )
                if ec_lc not in ['#888888', '#000000']:
                    print(seq.pep_seq, mod, ind + 1)
                ax.add_artist(
                    Ellipse(
                        (mod_col, mod_row - .04), 
                        1,
                        1,
                        ec=ec_lc,
                        fc='#ffffff',
                        lw=1.5,
                        zorder={
                            '#888888': 11,
                            '#000000': 11,
                        }.get(ec_lc, 12),
                    ),
                )
            
            if len(mods) > 0:
                fc_lc = '#cccccc'
            
            count[row, col] += 1

            # Draw peptide lines
            rarrow = True
            # print(seq, fc_lc)
            while col2 > 0:
                _draw_line(
                    ax,
                    row,
                    col, min([col2, max_col]),
                    rarrow=rarrow,
                    larrow=col2 < max_col,
                    line_color=fc_lc,
                    z={
                        '#cccccc': 1,
                        '#888888': 2,
                    }.get(fc_lc, 3),
                    row_offset=ro,
                )
                row += 1
                col = 0
                col2 -= max_col
                rarrow = False
        
    return figs