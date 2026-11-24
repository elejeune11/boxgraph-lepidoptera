import math
from scipy.optimize import brentq


def run(inputs, settings):
    kf_n = inputs["k_f_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]
    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0

    def f_bd(kb_try):
        p0 = k_f / (k_f + kb_try)
        return 16.0 * k_0 * (p0 ** 3) - kb_try

    kb_upper = None
    lo, hi = 1e-8 * k_0, 2.0 * k_0
    if f_bd(lo) * f_bd(hi) <= 0.0:
        kb_upper = brentq(f_bd, lo, hi, xtol=1e-13 * k_0, rtol=1e-11) / k_0

    def f_sp(kb_n_try):
        ksum = (kf_n + kb_n_try) * k_0
        rad = 1.0 - w * math.sqrt(ksum / D_B)
        hw = 0.0 if rad <= 0.0 else (ksum / math.sqrt(3.0)) * math.sqrt(rad)
        return hw - abs(kf_n - kb_n_try) * k_0

    sp_upper = None
    sp_lower = None
    if f_sp(kf_n) > 0.0:
        hi_u = max(kf_n * 1.0001, 1.0 - kf_n)
        if f_sp(hi_u) <= 0.0:
            sp_upper = brentq(f_sp, kf_n, hi_u, xtol=1e-12, rtol=1e-11)
        lo_l = 1e-9
        if lo_l < kf_n and f_sp(lo_l) <= 0.0:
            sp_lower = brentq(f_sp, lo_l, kf_n, xtol=1e-12, rtol=1e-11)

    return {
        "kb_Bdroplet_upper_norm": kb_upper,
        "kb_Bdroplet_lower_norm": 3.0 * kf_n,
        "kb_spinodal_upper_norm": sp_upper,
        "kb_spinodal_lower_norm": sp_lower,
    }
