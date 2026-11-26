"""Calibrated figure reader wrapper for fig4d_og_actin."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4bd

def read(image, at):
    return _fig4bd.read_series(image, 'og_actin', at)
