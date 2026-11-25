"""Calibrated reader for figs6b_surf_rs_0_1000: Fig. S6B Ring-Shell ratio 0:1000."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=3, row_idx=1, color_rgb=(184, 210, 210), color_tol=26,
        at=[float(v) for v in at], x_tick_vals=[0,100,200,300,400,500,600], y_tick_vals=[0.6,0.4,0.2,0.0],
        max_lum=218, min_sat=12, max_sat=48, max_res_px=5.0, label="Fig. S6B Ring-Shell 0:1000",
    )
