"""Fig. 8a calponin titration at 50 mM KCl."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig8 as M


def read(image, at):
    return M.read_fig8a(image, 'kcl_50mM')
