import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    R1_n = inputs["R1_norm"]
    kf_n = inputs["k_f_norm"]
    ratio = inputs["kb_over_kf"]
    Vs_n = inputs["Vs_norm"]
    phi_B_bar = inputs["phi_B_bar"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0
    k_b = ratio * k_f
    gb = w / 6.0
    V_s = Vs_n * (w ** 3)
    R1 = R1_n * w

    p0_req = gb / R1 + k_b * (R1 ** 2) / (3.0 * D_B)
    R2_norm = None
    R2_smallV_norm = None
    R_sym_norm = None
    R_sym_smallV_norm = None

    if kf_n == 0.0:
        V1 = (4.0 * math.pi / 3.0) * (R1 ** 3)
        if p0_req < 1.0:
            Vtot = V_s * (phi_B_bar - p0_req) / (1.0 - p0_req)
            V2 = Vtot - V1
            if V2 > 0.0:
                R2_norm = ((3.0 * V2 / (4.0 * math.pi)) ** (1.0 / 3.0)) / w
        Vtot_sv = V_s * (phi_B_bar - p0_req)
        V2_sv = Vtot_sv - V1
        if V2_sv > 0.0:
            R2_smallV_norm = ((3.0 * V2_sv / (4.0 * math.pi)) ** (1.0 / 3.0)) / w

        R_hi = (3.0 * V_s * phi_B_bar / (8.0 * math.pi)) ** (1.0 / 3.0) * 0.999999
        grid = np.geomspace(max(gb * 0.01, R_hi * 1e-6), R_hi, 300)
        f_ex = lambda r: ((V_s * phi_B_bar - 8.0 * math.pi / 3.0 * r ** 3) / (V_s - 8.0 * math.pi / 3.0 * r ** 3)) / r - gb / (r ** 2)
        vals = [f_ex(float(r)) for r in grid]
        rts = [brentq(f_ex, float(grid[i]), float(grid[i + 1]), xtol=1e-11 * w, rtol=1e-11) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0.0]
        if len(rts) >= 2:
            R_sym_norm = rts[-1] / w

        f_sv = lambda r: (phi_B_bar - 8.0 * math.pi / 3.0 * (r ** 3) / V_s) - gb / r
        r_pk = (3.0 * V_s * gb / (24.0 * math.pi)) ** 0.25
        if r_pk < R_hi and f_sv(r_pk) > 0.0:
            R_sym_smallV_norm = brentq(f_sv, r_pk, R_hi, xtol=1e-11 * w, rtol=1e-11) / w
    else:
        rhs = (k_f + 8.0 * math.pi * D_B * gb / V_s) / p0_req - (k_f + k_b)
        R2 = rhs * V_s / (4.0 * math.pi * D_B) - R1
        if R2 > 0.0:
            R2_norm = R2 / w
            R2_smallV_norm = R2_norm
        R_hi = max(math.sqrt(3.0 * D_B / k_b) * 2.0, (3.0 * V_s / (8.0 * math.pi)) ** (1.0 / 3.0))
        grid = np.geomspace(max(gb * 0.01, R_hi * 1e-6), R_hi, 300)
        f_ac = lambda r: D_B * (((k_f + 8.0 * math.pi * D_B * gb / V_s) / (k_f + k_b + 8.0 * math.pi * D_B * r / V_s)) / r - gb / (r ** 2)) - k_b * r / 3.0
        vals = [f_ac(float(r)) for r in grid]
        rts = [brentq(f_ac, float(grid[i]), float(grid[i + 1]), xtol=1e-11 * w, rtol=1e-11) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0.0]
        if len(rts) >= 2:
            R_sym_norm = rts[-1] / w
            R_sym_smallV_norm = R_sym_norm

    return {
        "R2_nullcline_norm": R2_norm,
        "R_sym_norm": R_sym_norm,
        "R2_nullcline_smallV_norm": R2_smallV_norm,
        "R_sym_smallV_norm": R_sym_smallV_norm,
    }
