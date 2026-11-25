"""Calibrated reader for fig4j_surf_kb1_ku0p01: Fig. 4J (k_bind=1.0, k_unbind=0.01 1/s)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ts_stacked as TS

def read(image, at):
    return TS.read_ts_stacked(
        image, n_rows=3, row_idx=1, color_rgb=(0, 0, 255), color_tol=55,
        at=[float(v) for v in at], x_tick_vals=[0, 100, 200, 300, 400, 500, 600],
        max_lum=130, min_sat=150, max_sat=None, label="Fig. 4J k_bind=1.0 k_unbind=0.01",
    )
