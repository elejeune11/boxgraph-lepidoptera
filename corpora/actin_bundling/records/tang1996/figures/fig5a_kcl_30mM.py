"""Fig. 5a Co(NH3)6^3+ titration at 30 mM KCl."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5 as M


def read(image, at):
    return M.read_fig5a(image, at, 'kcl_30mM')
