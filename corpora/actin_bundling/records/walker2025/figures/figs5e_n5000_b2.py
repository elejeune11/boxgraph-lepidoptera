"""Calibrated reader for figs5e_n5000_b2: Fig. S5E N=5000 dimers bound to 2 filaments."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=2, row_idx=0, x_vals=[0.001, 0.01, 0.1, 1.0, 10.0],
        palette_rgb=[(144, 204, 71), (200, 124, 35), (162, 47, 49)], cls_idx=2, y_tick_vals=[1.0, 0.5, 0.0],
        color_tol=60.0, max_res_px=5.0, label="Fig. S5E N=5000 bound=2",
    )
