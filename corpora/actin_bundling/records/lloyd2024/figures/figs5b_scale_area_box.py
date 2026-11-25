"""Digitize SI Fig. 5B: adult H. sara scale area boxplots (4 groups), PDF p. 28, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs5 import read_figs5_boxplots


def read(image, at):
    return read_figs5_boxplots(image, y_labels=[5000.0, 4000.0, 3000.0, 2000.0],
                               inset_row_max=960, inset_col_min=4480)
