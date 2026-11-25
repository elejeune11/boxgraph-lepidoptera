import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels as P

def read(image, at):
    return P.read_smfig5_flux(image, at, 5.0, -2.0, "diffusive", "d")
