import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    kf_n = inputs["k_f_norm"]
    ratio = inputs["kb_over_kf"]
    Vs_norm = inputs["Vs_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0
    k_b = ratio * k_f
    gb = w / 6.0
    V_s = Vs_norm * (w ** 3)

    def roots_m(m):
        def f(R):
            p0 = (k_f + 4.0 * math.pi * m * D_B * gb / V_s) / (k_f + k_b + 4.0 * math.pi * m * D_B * R / V_s)
            return D_B * (p0 / R - gb / (R ** 2)) - k_b * R / 3.0

        R_hi = max(math.sqrt(3.0 * D_B / k_b) * 2.0, (3.0 * V_s / (4.0 * math.pi * m)) ** (1.0 / 3.0))
        R_lo = max(gb * 0.01, R_hi * 1e-6)
        grid = np.geomspace(R_lo, R_hi, 300)
        vals = [f(float(r)) for r in grid]
        rts = []
        for i in range(len(grid) - 1):
            if vals[i] * vals[i + 1] <= 0.0:
                if vals[i] == 0.0:
                    rts.append(float(grid[i]))
                elif vals[i + 1] != 0.0:
                    rts.append(float(brentq(f, float(grid[i]), float(grid[i + 1]), xtol=1e-11 * R_lo, rtol=1e-10)))
        return sorted(rts)

    def is_stable_m(m):
        rm = roots_m(m)
        if len(rm) < 2:
            return False
        Rb = rm[-1]
        p0 = (k_f + 4.0 * math.pi * m * D_B * gb / V_s) / (k_f + k_b + 4.0 * math.pi * m * D_B * Rb / V_s)
        J00 = -(k_f + k_b) - 4.0 * math.pi * D_B * m * Rb / V_s
        Jii = D_B * (2.0 * gb / (Rb ** 3) - p0 / (Rb ** 2)) - k_b / 3.0
        if Jii >= 0.0:
            return False
        J0i = -4.0 * math.pi * D_B * p0 / V_s
        Ji0 = D_B / Rb
        # 2x2 collective block stability: tr < 0 and det > 0
        tr = J00 + Jii
        det = J00 * Jii - m * J0i * Ji0
        return (tr < 0.0) and (det > 0.0)

    r1 = roots_m(1)
    R_crit_norm = r1[0] / w if len(r1) >= 2 else None
    R_bar_norm = r1[-1] / w if len(r1) >= 2 else None

    # Binary search for maximal stable m in [1, 500]
    if not is_stable_m(1):
        N_max = 0
    else:
        lo, hi = 1, 2
        while hi < 500 and is_stable_m(hi):
            lo = hi
            hi = min(500, hi * 2)
        if hi == 500 and is_stable_m(500):
            N_max = 500
        else:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if is_stable_m(mid):
                    lo = mid
                else:
                    hi = mid
            N_max = lo

    N_max_chi1 = (V_s / (2.0 * math.pi * D_B * gb)) * (k_f * k_b / (k_f + k_b))
    return {
        "R_bar_norm": R_bar_norm,
        "R_crit_norm": R_crit_norm,
        "N_max": float(N_max),
        "N_max_chi1": N_max_chi1,
    }
