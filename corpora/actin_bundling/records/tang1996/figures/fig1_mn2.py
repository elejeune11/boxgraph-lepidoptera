"""Fig. 1 Mn2+ light scattering titration."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1 as M


def read(image, at):
    return M.read_fig1(image, at, 'mn2')
