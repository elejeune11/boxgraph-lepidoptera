"""Calibrated reader for fig5c_surf_kf1_ks0p1: Fig. 5C (k_dimer_form=1.0, k_dimer_split=0.1 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=2, row_idx=1, color_rgb=(102, 102, 102), color_tol=30,
        at=[float(v) for v in at], x_tick_vals=[0, 100, 200, 300, 400, 500, 600],
        max_lum=130, min_sat=None, max_sat=20, label="Fig. 5C k_form=1.0 k_split=0.1",
    )
