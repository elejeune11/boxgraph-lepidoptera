"""Figure reading wrapper delegating to _helpers.read_fig4_exp_by_shape."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _helpers as H


def read(image, at):
    return H.read_fig4_exp_by_shape(image, at, 'd225')
