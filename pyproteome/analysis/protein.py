from matplotlib.patches import Circle, Ellipse, Rectangle, FancyBboxPatch
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

def get_protein_seq(slc, gene):
    for _, row in slc.iterrows():
        for protein in row['Proteins'].proteins:
            if protein.gene == gene:
                return protein.full_sequence
            
def get_lc(fc, p, edge=False):
    if p < 1e-2:
        if fc > 1.25:
            return '#ff0000' if edge else '#ff8888'
        elif fc < 1/1.25:
            return '#0000ff' if edge else '#8888ff'
    return '#000000' if edge else '#888888'

def draw_line(
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
    lr_arrow = '{}-{}, head_width=.2, head_length=.2'.format(
        '<' if larrow else '',
        '>' if rarrow else '',
    )
    
    return ax.add_artist(
        FancyBboxPatch(
            (col + col_offset, row - .35), 
            col2 - col, 
            .6,
            ec=None,
            fc=line_color,
            zorder=-1,
            boxstyle="round,pad=.05",
        )
    )

    return ax.annotate(
        "",
        xy=(
            col + col_offset,
            row + row_offset,
        ),
        xytext=(
            col2 + col_offset,
            row + row_offset,
        ),
        xycoords='data',
        textcoords='data',
        arrowprops=dict(
            arrowstyle=lr_arrow,
            ec=line_color,
        ),
        zorder=z,
    )

def draw_protein_seq(
    ds,
    genes,
    max_col=50,
):
    ds = ds.filter(protein=genes)
    # ds = ds.filter(missed_cleavage=0)
    ds = ds.filter(missed_cleavage=1)

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

    print(mod_types, colors)

    for gene in genes:
        slc = ds[ds['Proteins'] == gene]
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
#         display(slc[['Proteins', 'Sequence', 'Fold Change', 'p-value']])

        z = 0
        prot_seq = get_protein_seq(slc, gene)
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
            )
            col += 1
            if col >= max_col:
                col = 0
                row += 1
                

        count = Counter()

        for _, pep in slc.iterrows():
            seq = pep['Sequence']
            match = [i for i in seq.protein_matches if i.protein.gene == gene][0]

            fc_lc = get_lc(pep['Fold Change'], pep['p-value'])
            ec_lc = get_lc(pep['Fold Change'], pep['p-value'], edge=True)

            mods = seq.modifications.get_mods(mod_types)
            print(seq.pep_seq, match.rel_pos, match.exact, str(seq))

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
                    color=ec_lc,
                    family='monospace',
                )
                ax.add_artist(
                    Ellipse(
                        (mod_col, mod_row - .04), 
                        1,
                        1,
                        ec=ec_lc,
                        fc='#ffffff',
                        lw=1.5,
#                         fc=colors[mod_types.index(mod.mod_type)],
                    ),
                )
            
            if len(mods) > 0:
                fc_lc = '#cccccc'
            
            count[row, col] += 1

            # Draw peptide lines
            rarrow = True
            while col2 > 0:
                draw_line(
                    ax,
                    row,
                    col, min([col2, max_col]),
                    rarrow=rarrow,
                    larrow=col2 < max_col,
                    line_color=fc_lc,
                    z=z,
                    row_offset=ro,
                )
                row += 1
                col = 0
                col2 -= max_col
                rarrow = False
                
            z += 1
        
    return figs