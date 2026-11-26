"""Fig. 4 Lys18 titration, no gelsolin."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4 as M


def read(image, at):
    return M.read_fig4(image, at, 'no_gel')
