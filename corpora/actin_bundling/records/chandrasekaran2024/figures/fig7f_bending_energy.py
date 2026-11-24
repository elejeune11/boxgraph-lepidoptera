import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig7f import read_fig7f
def read(image, at): return read_fig7f(image, "bending", at)
