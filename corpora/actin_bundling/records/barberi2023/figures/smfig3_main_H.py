import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels as P

def read(image, at):
    return P.read_smfig3_main(image, at, "H")
