"""Normalized cover scale bounding-box volume and surface area curves (% of max, Fig. 5F).

Evaluates V(t) = L(t)*W(t)*T(t) and S(t) = 2*(L*W + L*T + W*T) at the timepoints of L(t)
(interpolating W and T onto L's timepoints if needed) and normalizes each curve to 100%
at its maximum value.
"""

import numpy as np


def run(inputs, settings):
    cL, cW, cT = inputs["L"], inputs["W"], inputs["T"]
    tL = np.asarray(cL["x"], dtype=float)
    L = np.asarray(cL["y"], dtype=float)
    tW = np.asarray(cW["x"], dtype=float)
    W = np.asarray(cW["y"], dtype=float)
    tT = np.asarray(cT["x"], dtype=float)
    T = np.asarray(cT["y"], dtype=float)

    W_on_L = W if (tW.shape == tL.shape and np.allclose(tW, tL)) else np.interp(tL, tW, W)
    T_on_L = T if (tT.shape == tL.shape and np.allclose(tT, tL)) else np.interp(tL, tT, T)

    V = L * W_on_L * T_on_L
    S = 2.0 * (L * W_on_L + L * T_on_L + W_on_L * T_on_L)
    V_max = float(np.max(V)) if V.size else 1.0
    S_max = float(np.max(S)) if S.size else 1.0
    V_pct = 100.0 * V / V_max
    S_pct = 100.0 * S / S_max
    return {
        "V_pct": {"x": tL.tolist(), "y": V_pct.tolist()},
        "S_pct": {"x": tL.tolist(), "y": S_pct.tolist()},
    }
