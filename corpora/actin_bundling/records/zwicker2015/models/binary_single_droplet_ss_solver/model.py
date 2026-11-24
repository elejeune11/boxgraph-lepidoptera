import math
from scipy.optimize import brentq


def _peak_rate(D_B, gb, phi_0_B, k_b, rtol):
    # Peak of dR/dt = D_B*(phi_0_B/R - gb/R^2) - k_b*R/3 satisfies:
    # k_b*R^3 + 3*D_B*phi_0_B*R - 6*D_B*gb = 0
    if k_b <= 0.0:
        R_pk = 2.0 * gb / phi_0_B
    else:
        def dp(R):
            return k_b * R**3 + 3.0 * D_B * phi_0_B * R - 6.0 * D_B * gb

        R_hi = max(2.0 * gb / phi_0_B, (6.0 * D_B * gb / k_b) ** (1.0 / 3.0)) * 2.0
        R_pk = brentq(dp, 0.0, R_hi, xtol=rtol * gb * 1e-3, rtol=rtol)
    val_pk = D_B * (phi_0_B / R_pk - gb / R_pk**2) - k_b * R_pk / 3.0
    return R_pk, val_pk


def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    rtol = float(settings.get("rtol", 1.0e-6))

    phi_0_B = k_f / (k_f + k_b)
    ratio = k_b / k_f  # fixed phi_0_B = 1 / (1 + ratio)

    # Solve numerically for saddle-node k_f^c where peak dR/dt vanishes at fixed phi_0_B
    def peak_vs_kf(kf_try):
        kb_try = ratio * kf_try
        _, vpk = _peak_rate(D_B, gb, phi_0_B, kb_try, rtol * 0.1)
        return vpk

    kf_lo = 1e-6 * D_B / (gb**2) * phi_0_B**4
    while peak_vs_kf(kf_lo) <= 0.0:
        kf_lo *= 0.1
    kf_hi = max(kf_lo * 10.0, 10.0 * D_B / (gb**2) * phi_0_B**4 / max(1.0 - phi_0_B, 1e-9))
    while peak_vs_kf(kf_hi) >= 0.0:
        kf_hi *= 2.0
    k_f_c = brentq(peak_vs_kf, kf_lo, kf_hi, xtol=rtol * kf_lo, rtol=rtol)

    # Find exact roots R_crit and R_bar at given (k_f, k_b)
    R_pk, val_pk = _peak_rate(D_B, gb, phi_0_B, k_b, rtol * 0.1)
    if val_pk < 0.0:
        return {"R_crit": None, "R_bar": None, "k_f_c": k_f_c}

    def g(R):
        return D_B * (phi_0_B / R - gb / R**2) - k_b * R / 3.0

    R_lo = gb / phi_0_B * 0.1
    while g(R_lo) >= 0.0:
        R_lo *= 0.5
    R_crit = brentq(g, R_lo, R_pk, xtol=rtol * R_lo, rtol=rtol)

    R_hi = max(R_pk * 2.0, math.sqrt(3.0 * D_B * phi_0_B / k_b) * 2.0)
    while g(R_hi) >= 0.0:
        R_hi *= 2.0
    R_bar = brentq(g, R_pk, R_hi, xtol=rtol * R_pk, rtol=rtol)

    return {"R_crit": R_crit, "R_bar": R_bar, "k_f_c": k_f_c}
