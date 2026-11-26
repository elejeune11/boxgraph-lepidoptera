"""Fig. 7 MARCKS peptide titration curve."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig7 as M


def read(image, at):
    return M.read_fig7(image, 'marcks')
