"""Fig. 1 His5 light scattering titration."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1 as M


def read(image, at):
    return M.read_fig1(image, at, 'his5')
