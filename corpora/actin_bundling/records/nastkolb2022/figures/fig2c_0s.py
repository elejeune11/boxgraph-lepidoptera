"""Digitize Fig. 2C transverse actin bundle fluorescence intensity profile at 0 s (darkest purple; PDF p. 4, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2c_common import read_fig2c

def read(image, at):
    return read_fig2c(image, "0s", at)
