"""Calibrated reader for figs7d_surf_kf0p001_ks10: Fig. S7D (k_multimer_form=0.001, physical k_multimer_split=10.0 1/s, pale teal (mislabeled 0.001 in inverted legend))."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=5, row_idx=4, color_rgb=(184, 210, 210), color_tol=26,
        at=[float(v) for v in at], x_tick_vals=[0,100,200,300,400,500,600], y_tick_vals=[0.6,0.3,0.0],
        max_lum=218, min_sat=12, max_sat=48, max_res_px=5.0, label="Fig. S7D k_form=0.001 phys_k_split=10.0",
    )
