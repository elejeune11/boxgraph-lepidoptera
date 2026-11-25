"""Calibrated reader for figs1c_partition_minilpd: Fig. S1C mini-Lpd partition coefficient."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bar_replicates as BR

def read(image, at):
    return BR.read_bar_tops_and_replicates(image, [60,50,40,30,20,10,0], [50.0, 150.0], n_bars=4, select_indices=[0, 2], label="Fig. S1C mini-Lpd")
