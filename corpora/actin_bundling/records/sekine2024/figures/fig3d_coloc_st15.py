"""Digitize Stage 15 Pearson's R colocalization boxplots (7 proteins: GFP, DAAM, alpha-Actinin, Zasp52, Myosin II, ROCK, RhoGEF2) in Fig. 3d (PDF p. 5, 900 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _coloc_boxplot_common import read_coloc_boxplots


def read(image, at):
    return read_coloc_boxplots(image, n_expected=7, min_bar_w=100, min_bar_h=8)
