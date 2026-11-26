"""Fig. 2b Mg2+ titration with no Cu2+."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2 as M


def read(image, at):
    return M.read_fig2b(image, at, 'mg_no_cu')
