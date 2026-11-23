"""Standard 2D shape circularity C = 4 * pi * A / P^2 (README.md line 5; sekine2024 Methods p. 10).

Note: Nanocluster_motion.m never computes perimeter or circularity. Also provides discrete binary
mask circularity helpers comparing outer-edge pixel perimeter (4-edge boundary length, matching
rectilinear ImageJ pixel boundary on axis-aligned shapes) and MATLAB regionprops('Circularity')
small-region correction 4*pi*A / (P_center + pi)^2.
"""

import math
import numpy as np


def edge_perimeter(mask):
    """Outer pixel-edge perimeter (count of foreground-background 4-edges) of a binary mask."""
    m = np.asarray(mask, dtype=np.int32)
    p = np.pad(m, 1, mode="constant", constant_values=0)
    return float(np.sum(np.abs(np.diff(p, axis=0))) + np.sum(np.abs(np.diff(p, axis=1))))


def circularity_from_mask(mask):
    """Circularity 4*pi*A/P^2 using outer pixel-edge perimeter (clamped at 1.0 as in ImageJ)."""
    area = float(np.asarray(mask, dtype=bool).sum())
    perim = edge_perimeter(mask)
    if area <= 0.0 or perim <= 0.0:
        return 0.0
    return min(1.0, 4.0 * math.pi * area / (perim * perim))


def run(inputs, settings):
    A = float(inputs["A"])
    P = float(inputs["P"])
    return {"C": 4.0 * math.pi * A / (P * P)}


def symbolic():
    import sympy as sp

    A, P = sp.symbols("A P", positive=True)
    return {"C": 4 * sp.pi * A / P**2}
