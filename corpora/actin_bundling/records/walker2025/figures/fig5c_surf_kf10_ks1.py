"""Calibrated reader for fig5c_surf_kf10_ks1: Fig. 5C (k_dimer_form=10, k_dimer_split=1.0 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=2, row_idx=0, color_rgb=(198, 198, 198), color_tol=18,
        at=[float(v) for v in at], x_tick_vals=[0, 100, 200, 300, 400, 500, 600],
        max_lum=210, min_sat=None, max_sat=15, label="Fig. 5C k_form=10 k_split=1.0",
    )
