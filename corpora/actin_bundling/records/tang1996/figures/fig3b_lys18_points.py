"""Fig. 3b 6 measured Lys18 onset points vs [Actin]."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3 as M


def read(image, at):
    return M.read_fig3b_points(image)
