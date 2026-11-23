"""Compute the exact interval (dx_lo, dx_hi] of square pixel edge lengths dx for which
[N_min, N_max] (Nanocluster_motion.m line 40: [6, 53] px) is the exact set of integer pixel
areas inside [A_min, A_max] (README.md line 5: [0.01, 0.1] um^2):
  (N_min - 1) * dx^2 < A_min <= N_min * dx^2
  N_max * dx^2 <= A_max < (N_max + 1) * dx^2
"""

import math


def run(inputs, settings):
    n_min = float(inputs["N_min"])
    n_max = float(inputs["N_max"])
    a_min = float(inputs["A_min"])
    a_max = float(inputs["A_max"])
    dx_lo = max(math.sqrt(a_min / n_min), math.sqrt(a_max / (n_max + 1.0)))
    dx_hi = min(math.sqrt(a_min / (n_min - 1.0)), math.sqrt(a_max / n_max))
    return {
        "dx_lo": dx_lo,
        "dx_hi": dx_hi,
        "dx_mid": 0.5 * (dx_lo + dx_hi),
    }
