import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5d import read_fig5d
def read(image, at): return read_fig5d(image, 0)
