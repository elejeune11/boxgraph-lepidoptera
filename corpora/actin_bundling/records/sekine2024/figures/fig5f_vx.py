"""Digitize Early stage 16 PIV longitudinal mean displacement Vx boxplots (1=Control, 2=DAAM RNAi, 3=Zasp52 RNAi, 4=myosin II DN) in Fig. 5f (PDF p. 9, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5f_common import read_fig5f


def read(image, at):
    return read_fig5f(image, "vx")
