"""Calibrated figure reader for figS1a_B."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1_S3 import read_figS1

def read(image, at):
    return read_figS1(image, at, "a", "B")
