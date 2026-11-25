"""Calibrated reader for fig5d_kf1_singly_bound_dimer: Fig. 5D (k_dimer_form=1.0 1/s, Singly Actin-Bound Dimers)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=2, row_idx=1, x_vals=[10.0, 1.0, 0.1, 0.01, 0.001],
        palette_rgb=[(1, 0, 128), (128, 0, 0), (103, 103, 103), (182, 209, 209), (9, 121, 105)], cls_idx=3, y_tick_vals=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
        label="Fig. 5D k_form=1.0 Singly Actin-Bound Dimers",
    )
