import numpy as np
from scipy.optimize import brentq


def _max_growth_rate(d, rho_a, rho_s, Pe, qs):
    q2 = qs * qs
    act = 0.25 * Pe * q2 / (1.0 + q2)
    a11 = -d * q2 + rho_a + act
    a12 = rho_a
    a21 = -2.0 * rho_s + act
    a22 = -q2 - rho_s
    tr = a11 + a22
    det = a11 * a22 - a12 * a21
    disc = tr * tr - 4.0 * det
    re = np.where(disc >= 0.0, 0.5 * (tr + np.sqrt(np.maximum(disc, 0.0))), 0.5 * tr)
    idx = int(np.argmax(re))
    q_best = float(qs[idx])
    re_best = float(re[idx])
    # Parabolic interpolation around interior grid maximum for smooth O(dk^4) accuracy
    if 0 < idx < len(qs) - 1:
        y0, y1, y2 = float(re[idx - 1]), float(re[idx]), float(re[idx + 1])
        denom = y0 - 2.0 * y1 + y2
        if denom < -1e-18:
            delta = 0.5 * (y0 - y2) / denom
            if abs(delta) <= 1.0:
                dq = float(qs[1] - qs[0])
                q_best = float(qs[idx] + delta * dq)
                re_best = float(y1 - 0.25 * (y0 - y2) * delta)
    return re_best, q_best


def run(inputs, settings):
    rho_a = float(inputs["rho_a"])
    rho_s = float(inputs["rho_s"])
    Pe = float(inputs["Pe"])
    n_k = int(round(settings.get("n_k", 200)))
    k_max = float(settings.get("k_max_search", 20.0))

    if rho_s <= rho_a:
        raise ValueError("Region I/II Turing boundary requires rho_s > rho_a (homogeneous stability at k = 0)")

    # Exclude k = 0 (strictly negative when rho_s > rho_a)
    qs = np.linspace(k_max / n_k, k_max, n_k)

    def f_d(d):
        val, _ = _max_growth_rate(d, rho_a, rho_s, Pe, qs)
        return val

    d_lo = 1e-6
    d_hi = max(2.0, 0.5 * Pe + 2.0)
    while f_d(d_hi) > 0.0:
        d_hi *= 2.0
        if d_hi > 1e6:
            raise RuntimeError("Could not bracket critical Da/Ds")

    d_crit = float(brentq(f_d, d_lo, d_hi, xtol=1e-12, rtol=1e-12))
    _, k_crit = _max_growth_rate(d_crit, rho_a, rho_s, Pe, qs)
    return {"Da_over_Ds_crit": d_crit, "k_crit_l": float(k_crit)}
