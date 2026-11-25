"""Calibrated reader for figs6g_dimers_kb0p001_b1: Fig. S6G (k_bind=0.001 1/s, bound=1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=5, row_idx=4, x_vals=[0.001, 0.01, 0.1, 1.0, 10.0],
        palette_rgb=[(144, 204, 71), (200, 124, 35), (162, 47, 49)], cls_idx=1, y_tick_vals=[1.0, 0.5, 0.0],
        color_tol=60.0, max_res_px=5.0, label="Fig. S6G k_bind=0.001 bound=1",
    )
