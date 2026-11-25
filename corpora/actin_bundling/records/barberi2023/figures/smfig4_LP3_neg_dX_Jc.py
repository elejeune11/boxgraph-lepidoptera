import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels as P

def read(image, at):
    return P.read_smfig4(image, at, "neg_dX_Jc", "LP3")
