"""Fig. 6b ATP bundle dissociation curve."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig6 as M


def read(image, at):
    return M.read_fig6b(image, at, 'atp')
