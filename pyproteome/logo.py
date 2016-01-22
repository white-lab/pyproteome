
# Built-ins
from math import ceil

# IPython
from IPython.display import display

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np


def draw_logo(n_mers):
    """
    Draw a sequence logo.

    Parameters
    ----------
    n_mers : list of str
    """
    assert len(n_mers) >= 1

    n = len(n_mers[0])
    assert n % 2 == 1
    assert all(len(i) == n for i in n_mers)

    f, ax = plt.subplots()

    indices = np.arange(-n // 2, n // 2 + 1)
    bits = [{"A": 1, "S": 2}]

    ax.set_xlabel("Position")
    ax.set_xticks(indices)
    #ax.set_xticklabels

    ax.set_ylabel("Bits")
    max_bits = int(ceil(max(max(i.values()) for i in bits)))
    ax.set_yticks(np.arange(max_bits + 1))
    ax.set_ylim((0, max_bits + .5))

    from matplotlib.transforms import Affine2D

    for index, seq_bits in enumerate(bits):
        x_pos = indices[index]
        y_pos = 0

        for letter, size in seq_bits.items():
            tr = Affine2D().scale(1, 2)
            display(letter)

            plt.text(x_pos, y_pos, letter, transform=ax.transData)
            y_pos += size

    return f
