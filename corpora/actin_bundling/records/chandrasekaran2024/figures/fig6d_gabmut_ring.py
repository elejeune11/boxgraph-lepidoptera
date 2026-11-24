import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig6d import read_fig6d
def read(image, at): return read_fig6d(image, "GABmut", "ring")
