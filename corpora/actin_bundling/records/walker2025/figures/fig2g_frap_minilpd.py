import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2g import read_fig2g

def read(image, at):
    return read_fig2g(image, "minilpd", at)
