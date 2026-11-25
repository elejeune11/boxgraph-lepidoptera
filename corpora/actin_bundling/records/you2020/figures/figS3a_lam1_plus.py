"""Calibrated figure reader for figS3a_lam1_plus."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1_S3 import read_figS3

def read(image, at):
    return read_figS3(image, at, "a", "plus")
