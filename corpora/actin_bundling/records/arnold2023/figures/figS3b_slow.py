"""Fig. S3(B) 'Slow' region (pink): a/a_0 vs. shifted time t [s]."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS3b import read_series


def read(image, at):
    return read_series(image, at, "slow")
