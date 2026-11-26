"""Figure reading wrapper delegating to _helpers.read_figS2a_series."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _helpers as H


def read(image, at):
    return H.read_figS2a_series(image, at, 32)
