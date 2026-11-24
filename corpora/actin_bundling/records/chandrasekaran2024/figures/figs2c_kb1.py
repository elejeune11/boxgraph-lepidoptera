import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2c import read_figs2c
def read(image, at): return read_figs2c(image, 1, "1.0/s")
