import math
from scipy.optimize import brentq


def run(inputs, settings):
    kf_n = inputs["k_f_norm"]
    ratio = inputs["kb_over_kf"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0
    k_b = ratio * k_f
    gb = w / 6.0
    phi_0_B = 1.0 / (1.0 + ratio)

    # Eq. (16), approx R_crit, and Eq. (17)
    r16 = math.sqrt(3.0 * D_B * phi_0_B / k_b)
    rc_app = gb / phi_0_B
    kbc = (4.0 / 9.0) * (D_B / (gb ** 2)) * (phi_0_B ** 3)
    kfc = kbc * phi_0_B / (1.0 - phi_0_B)
    R_star = (3.0 * D_B * gb / (2.0 * kbc)) ** (1.0 / 3.0)

    # Exact cubic roots of k_b*R^3 - 3*D_B*phi_0_B*R + 3*D_B*gb = 0
    dp = lambda R: k_b * (R ** 3) + 3.0 * D_B * phi_0_B * R - 6.0 * D_B * gb
    R_pk_hi = max(2.0 * gb / phi_0_B, (6.0 * D_B * gb / k_b) ** (1.0 / 3.0)) * 2.0
    R_pk = brentq(dp, 0.0, R_pk_hi, xtol=1e-13 * gb, rtol=1e-11)
    g = lambda R: D_B * (phi_0_B / R - gb / (R ** 2)) - k_b * R / 3.0
    if g(R_pk) < 0.0:
        R_crit = None
        R_bar = None
    else:
        R_lo = 0.1 * gb / phi_0_B
        while g(R_lo) >= 0.0:
            R_lo *= 0.5
        R_crit = brentq(g, R_lo, R_pk, xtol=1e-12 * R_lo, rtol=1e-11)
        R_hi = max(2.0 * R_pk, 2.0 * r16)
        while g(R_hi) >= 0.0:
            R_hi *= 2.0
        R_bar = brentq(g, R_pk, R_hi, xtol=1e-12 * R_pk, rtol=1e-11)

    return {
        "R_bar_norm": None if R_bar is None else R_bar / w,
        "R_crit_norm": None if R_crit is None else R_crit / w,
        "R_bar_eq16_norm": r16 / w,
        "R_crit_approx_norm": rc_app / w,
        "R_star_norm": R_star / w,
        "k_f_c_norm": kfc / k_0,
    }
