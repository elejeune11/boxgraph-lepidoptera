"""Calibrated figure reader for figS2_im_plus."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1_S3 import read_figS2

def read(image, at):
    return read_figS2(image, at, "im_plus")
