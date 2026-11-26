"""Calibrated reader for Fig. 6B (p. 9): -VASP (x 2.4) dashed reference curve."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6b(image)
    m = C.gray(image) < 160.0
    # Exclude upper-right legend and top-left '[KCl]=50mM' text
    m[: int(ya.to_pixel(4800.0)), int(xa.to_pixel(6.0)) :] = False
    m[: int(ya.to_pixel(6000.0)), :] = False
    # At x=1, 2, 4, the dashed curve lies strictly between the bottom -VASP curve and upper +50nM VASP curve
    def pick(x, cl):
        # Filter to clusters strictly below the +50nM VASP curve (top) and above -VASP (bottom)
        if len(cl) < 3:
            return None
        # At x=1 the hollow open circle splits into 2 top edges; dashed line is immediately above bottom -VASP cluster
        # Specifically, among clusters with y < 3000, the dashed curve is the second from bottom (-2)
        if x <= 1.5:
            # Find cluster immediately above bottom (-VASP ~619), around 1400-1800
            below_3k = [i for i, (rc, _, _) in enumerate(cl) if ya.to_value(rc) < 2800.0]
            if len(below_3k) >= 3:
                return below_3k[-3] if ya.to_value(cl[below_3k[-2]][0]) < 1200.0 else below_3k[-2]
            return below_3k[0] if len(below_3k) >= 2 else None
        # At x=2 and x=4, dashed curve is cleanly second from bottom (-2)
        return -2
    return C.read_ranked_dark_curves(
        image, xa, ya, fr, at, pick,
        mask=m, half_w=10, merge_gap=6, y_min_val=50.0, y_max_val=5800.0,
    )
