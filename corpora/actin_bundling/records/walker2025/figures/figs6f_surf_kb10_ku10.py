"""Calibrated reader for figs6f_surf_kb10_ku10: Fig. S6F (k_bind=10.0, k_unbind=10.0 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=5, row_idx=0, color_rgb=(184, 210, 210), color_tol=26,
        at=[float(v) for v in at], x_tick_vals=[0,100,200,300,400,500,600], y_tick_vals=[0.6,0.3,0.0],
        max_lum=218, min_sat=12, max_sat=48, max_res_px=5.0, label="Fig. S6F k_bind=10.0 k_unbind=10.0",
    )
