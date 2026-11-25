"""Calibrated reader for figs6b_surf_rs_250_750: Fig. S6B Ring-Shell ratio 250:750."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=3, row_idx=1, color_rgb=(198, 198, 198), color_tol=20,
        at=[float(v) for v in at], x_tick_vals=[0,100,200,300,400,500,600], y_tick_vals=[0.6,0.4,0.2,0.0],
        max_lum=212, min_sat=None, max_sat=15, max_res_px=5.0, label="Fig. S6B Ring-Shell 250:750",
    )
