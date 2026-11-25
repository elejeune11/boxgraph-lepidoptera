"""Calibrated reader for fig4k_dimers_kb1_b2: Fig. 4K (k_bind=1.0 1/s, 2 bound filaments)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=3, row_idx=1, x_vals=[0.01, 0.1, 1.0],
        palette_rgb=[(144, 204, 71), (200, 124, 35), (162, 47, 49)], cls_idx=2, y_tick_vals=[1.0, 0.5, 0.0],
        label="Fig. 4K k_bind=1.0 bound=2",
    )
