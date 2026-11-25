"""Calibrated reader for fig4j_surf_kb10_ku1: Fig. 4J (k_bind=10, k_unbind=1.0 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=3, row_idx=0, color_rgb=(102, 102, 102), color_tol=35,
        at=[float(v) for v in at], x_tick_vals=[0, 100, 200, 300, 400, 500, 600],
        max_lum=130, min_sat=None, max_sat=20, label="Fig. 4J k_bind=10 k_unbind=1.0",
    )
