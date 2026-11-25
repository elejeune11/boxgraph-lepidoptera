"""Calibrated reader for figs6d_lpd_rs_b0: Fig. S6D Ring-Shell mini-Lpd dimers bound to 0 filaments."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=3, row_idx=1, x_vals=[250.0, 500.0, 750.0, 1000.0],
        palette_rgb=[(144, 204, 71), (200, 124, 35), (162, 47, 49)], cls_idx=0, y_tick_vals=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
        active_bar_indices=[1, 2, 3, 4], n_slots=5, color_tol=60.0, max_res_px=5.0,
        label="Fig. S6D Ring-Shell mini-Lpd bound=0",
    )
