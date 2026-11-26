"""Fig. 3b plotted linear fit line [Lys18] = 19.4 + 0.26*[Actin]."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3 as M


def read(image, at):
    return M.read_fig3b_fit_line(image, at)
