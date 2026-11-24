import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs11d import read_figs11d
def read(image, at): return read_figs11d(image, 0, "ring_shell_sb")
