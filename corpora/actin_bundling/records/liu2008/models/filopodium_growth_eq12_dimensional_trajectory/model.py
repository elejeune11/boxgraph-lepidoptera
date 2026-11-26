"""Evaluate dimensionally restored SI Eq. (12) along target timestamps using measured c_0(t) from Fig. s5."""

from __future__ import annotations

import math
import numpy as np


def _integrate_c0(c0_x: np.ndarray, c0_y: np.ndarray, t0: float, t1: float, n_substeps: int) -> float:
    if t1 <= t0:
        return 0.0
    n = max(int(n_substeps), 2)
    grid = np.linspace(t0, t1, n + 1)
    vals = np.interp(grid, c0_x, c0_y, left=float(c0_y[0]), right=float(c0_y[-1]))
    return float(np.trapezoid(vals, grid))


def run(inputs, settings):
    c0 = inputs["c0_curve"]
    tgt = inputs["target_curve"]
    c0_x = np.asarray(c0["x"], dtype=float)
    c0_y = np.asarray(c0["y"], dtype=float)
    order = np.argsort(c0_x)
    c0_x, c0_y = c0_x[order], c0_y[order]

    tx = [float(v) for v in tgt["x"]]
    ty = [float(v) for v in tgt["y"]]
    if not tx:
        return {"L_mem_curve": {"x": [], "y": []}}

    t0 = tx[0]
    L0 = ty[0]
    f_mem = float(inputs["f_mem"])
    delta = float(inputs["delta"])
    N = float(inputs["N"])
    kBT = float(inputs["kBT"])
    L_star = float(inputs["L_star"])
    v_0 = float(inputs["v_0"])
    n_sub = int(settings.get("n_substeps", 100))

    bolt = math.exp(-f_mem * delta / (N * kBT))
    L0_dim = L0 / L_star
    cum_int = 0.0
    out_y = []
    prev_t = t0
    for t in tx:
        if t > prev_t:
            cum_int += _integrate_c0(c0_x, c0_y, prev_t, t, n_sub)
            prev_t = t
        rad = 1.0 + L0_dim * L0_dim + 2.0 * L0_dim + (2.0 * v_0 / L_star) * bolt * cum_int
        out_y.append(L_star * (math.sqrt(rad) - 1.0))
    return {"L_mem_curve": {"x": tx, "y": out_y}}
