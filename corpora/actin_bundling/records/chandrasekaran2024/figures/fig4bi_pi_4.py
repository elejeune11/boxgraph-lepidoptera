import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4b import read_fig4b
def read(image, at): return read_fig4b(image, "i", 0, "pi/4", at)
