"""Calibrated reader for figs7f_multimer_len_kf10: Fig. S7F median multimer length with actin (k_multimer_form=10.0 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _violin_si as VS

def read(image, at):
    return VS.read_violin_panel(image, n_rows=5, row_idx=0, mode="quartile_dash", label="Fig. S7F k_form=10.0")
