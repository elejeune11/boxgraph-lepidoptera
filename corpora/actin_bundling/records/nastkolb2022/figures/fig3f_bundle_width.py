"""Digitize Fig. 3F bundle width boxplots across 5 protein conditions (1=VASP, 2=VASP_DeltaTetra, 3=alpha-actinin, 4=mDia1, 5=mDia1+alpha-actinin; PDF p. 5, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3f_common import read_fig3f

def read(image, at):
    return read_fig3f(image, "width")
