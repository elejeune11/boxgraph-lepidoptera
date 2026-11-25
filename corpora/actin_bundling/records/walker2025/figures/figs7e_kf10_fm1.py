"""Calibrated reader for figs7e_kf10_fm1: Fig. S7E (k_multimer_form=10.0 1/s, state FM1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_stacked as BS

def read(image, at):
    return BS.read_stacked_bars_panel(
        image, n_rows=5, row_idx=0, x_vals=[10.0, 1.0, 0.1, 0.01, 0.001],
        palette_rgb=[(6, 6, 130), (126, 6, 6), (102, 102, 102), (210, 126, 46), (186, 210, 210), (10, 122, 102)], cls_idx=2, color_tol=55.0, max_res_px=5.0,
        use_box_y_fallback=True, label="Fig. S7E k_form=10.0 FM1",
    )
