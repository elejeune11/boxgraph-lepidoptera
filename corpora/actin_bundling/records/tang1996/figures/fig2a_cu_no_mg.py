"""Fig. 2a Cu2+ titration without MgCl2."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2 as M


def read(image, at):
    return M.read_fig2a(image, 'cu_no_mg')
