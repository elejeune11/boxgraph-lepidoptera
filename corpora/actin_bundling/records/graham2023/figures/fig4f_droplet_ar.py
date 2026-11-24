"""Calibrated geometric reader for Main Fig. 4f (22 simulated droplets at D = 1 um, gamma_I = 5.0 pN).

Measures the horizontal major span (w) and vertical minor span (h) of each of the 22 rendered
droplet ellipses (T_ring = 5 to 480 nm) directly from the rendered image. The ratio a/b = w/h
is scale-invariant; its uncertainty reflects the finite stroke width of the green actin ring
and the slight left cusp protrusion (~4 px).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

# Printed T_ring labels above (top row) and below (bottom row) the 22 subpanels in Fig. 4f (nm)
T_RING_LABELS_NM = [
    5.0, 30.0, 45.0, 75.0, 90.0, 120.0, 135.0, 165.0, 180.0, 210.0, 225.0,
    255.0, 270.0, 300.0, 315.0, 345.0, 360.0, 390.0, 405.0, 435.0, 450.0, 480.0,
]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    # Detect colored droplet+ring regions (green ring + teal droplet interior)
    teal_green = (g - r > 10.0) & (g > 60.0)
    lab, n = ndimage.label(ndimage.binary_dilation(teal_green, iterations=4))
    objs = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) > 4000:
            objs.append((float(rr.mean()), float(cc.mean()), int(rr.min()), int(rr.max()), int(cc.min()), int(cc.max())))
    if len(objs) != 22:
        raise RuntimeError(f"Expected 22 simulated droplets in Fig. 4f, found {len(objs)}")

    # Sort into top row (11) then bottom row (11), left to right
    y_split = 0.5 * (min(o[0] for o in objs) + max(o[0] for o in objs))
    objs.sort(key=lambda o: (0 if o[0] < y_split else 1, o[1]))

    H, W = image.shape[:2]
    all_readings = {}
    for tr_nm, o in zip(T_RING_LABELS_NM, objs):
        r0 = max(0, o[2] - 18)
        r1 = min(H, o[3] + 18)
        c0 = max(0, o[4] - 18)
        c1 = min(W, o[5] + 18)
        patch = image[r0:r1, c0:c1]
        # Full droplet extent (black ellipse boundary + teal fill + green ring)
        nw = patch.mean(axis=2) < 235.0
        pl, pn = ndimage.label(nw)
        sizes = [int(np.sum(pl == k)) for k in range(1, pn + 1)]
        kmax = int(np.argmax(sizes)) + 1
        pr, pc = np.nonzero(pl == kmax)
        w = float(pc.max() - pc.min() + 1)
        h = float(pr.max() - pr.min() + 1)
        ar = w / h
        # Conservative ±5 px uncertainty on w and ±2 px on h for left cusp & stroke thickness
        ar_lo = (w - 5.0) / (h + 2.0)
        ar_hi = (w + 2.0) / (h - 2.0)
        cx = c0 + 0.5 * (pc.min() + pc.max()) + 0.5
        cy = r0 + 0.5 * (pr.min() + pr.max()) + 0.5
        all_readings[tr_nm] = {
            "x": tr_nm,
            "y": float(ar),
            "lo": float(min(ar_lo, ar)),
            "hi": float(max(ar_hi, ar)),
            "px": [float(cx), float(cy)],
            "px_lo": float(r0 + pr.min()),
            "px_hi": float(r0 + pr.max()),
            "note": f"Fig4f droplet bbox w={w:.0f}px, h={h:.0f}px, sqrt(wh)={np.sqrt(w*h):.1f}px",
        }

    targets = at if at else T_RING_LABELS_NM
    out = []
    for x in targets:
        nearest = min(T_RING_LABELS_NM, key=lambda v: abs(v - x))
        if abs(nearest - x) <= 1.0:
            out.append(all_readings[nearest])
        else:
            out.append({"x": float(x), "y": None, "found": False, "note": f"no Fig. 4f subpanel at T_ring={x}"})
    return out
