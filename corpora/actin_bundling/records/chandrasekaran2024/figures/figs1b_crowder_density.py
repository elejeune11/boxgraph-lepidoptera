import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs1b import read_figs1b
def read(image, at): return read_figs1b(image, "crowder")
