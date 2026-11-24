import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig6b import read_fig6b
def read(image, at): return read_fig6b(image, "shell")
