"""Digitize mean longitudinal (0-30 deg) fusion duration across 4 genotypes in Fig. 3e (1=Control, 2=DAAM RNAi, 3=zasp52 RNAi, 4=myosin II DN, PDF p. 5, 900 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fusion_boxplot_common import read_fusion_means


def read(image, at):
    return read_fusion_means(image, [1.0, 2.0, 3.0, 4.0], angle_slot=2)
