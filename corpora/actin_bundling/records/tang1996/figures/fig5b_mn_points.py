"""Fig. 5b 4 measured Mn2+ onset points vs [KCl]."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5 as M


def read(image, at):
    return M.read_fig5b_mn_points(image)
