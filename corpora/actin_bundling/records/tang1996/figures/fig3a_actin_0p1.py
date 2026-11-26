"""Fig. 3a Lys18 titration at 0.1 mg/ml F-actin (3 replicates)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3 as M


def read(image, at):
    return M.read_fig3a(image, at, 'actin_0p1')
