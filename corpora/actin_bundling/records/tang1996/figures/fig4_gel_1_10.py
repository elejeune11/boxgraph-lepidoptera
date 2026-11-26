"""Fig. 4 Lys18 titration, 1:10 gelsolin:actin."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4 as M


def read(image, at):
    return M.read_fig4(image, at, 'gel_1_10')
