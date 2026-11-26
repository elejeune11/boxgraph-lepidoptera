"""Calibrated reader for Fig. S2C top (p. 26): normalized Actin line scan (actin alone)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS2c(image, at, "top", "actin")
