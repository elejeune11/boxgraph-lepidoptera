"""Digitize SI Fig. 5A: adult H. sara ridge spacing boxplots (4 groups), PDF p. 28, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs5 import read_figs5_boxplots


def read(image, at):
    return read_figs5_boxplots(image, y_labels=[1.4, 1.2, 1.0, 0.8, 0.6],
                               inset_row_max=700, inset_col_min=4480)
