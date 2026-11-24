import sys
from pathlib import Path
_d = str(Path(__file__).resolve().parent)
if _d not in sys.path:
    sys.path.insert(0, _d)
from _fig14_16 import read_fig15b_profile

def read(image, at):
    return read_fig15b_profile(image, at, "phi2")
