"""Fig. 5b Co(NH3)6^3+ power-law fit line."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5 as M


def read(image, at):
    return M.read_fig5b_fit_line(image, at, 'cohex')
