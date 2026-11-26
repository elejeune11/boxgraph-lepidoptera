import cmath
import math
from scipy.optimize import minimize_scalar

SITES_28_13 = ((0, 0), (4, 1), (9, -1), (14, 3), (19, 1), (24, -1))


def u_norm_periodic(xi_t: float, w0_signed: float) -> float:
    """Exact infinite-chain periodic geometric sum for u / (beta*U)^2."""
    z_p = cmath.exp(-28.0 / xi_t + 1j * 28.0 * w0_signed)
    total = 0.0
    for l1, m1 in SITES_28_13:
        for l2, m2 in SITES_28_13:
            dphi = 2.0 * math.pi * (m1 - m2) / 6.0
            dl0 = l1 - l2
            d_pos = dl0 % 28
            if d_pos == 0:
                d_pos = 28
            s_pos = (cmath.exp(-d_pos / xi_t + 1j * (w0_signed * d_pos - dphi)) / (1.0 - z_p)).real
            d_neg = (-dl0) % 28
            if d_neg == 0:
                d_neg = 28
            s_neg = (cmath.exp(-d_neg / xi_t + 1j * (w0_signed * d_neg + dphi)) / (1.0 - z_p)).real
            total += 0.5 * (s_pos + s_neg)
    return total / float(len(SITES_28_13))


def run(inputs, settings):
    beta_U = float(inputs["beta_U"])
    xi_t = float(inputs["xi_t"])
    # Apply left-handed sign (-omega_0) matching the chirality of SITES_28_13
    w0_lh = -abs(float(inputs["omega_0"]))
    opt_tol = float(settings.get("opt_tol", 1.0e-6))

    u_n = u_norm_periodic(xi_t, w0_lh)
    u_val = (beta_U ** 2) * u_n

    # Log-slope with respect to beta_U
    b1, b2 = beta_U * 0.99, beta_U * 1.01
    stiffness_exp = math.log(((b2 ** 2) * u_n) / ((b1 ** 2) * u_n)) / math.log(b2 / b1)

    # Local log-slope d(ln u)/d(ln xi_t) at input xi_t
    eps = 1.0e-4 * xi_t
    u_lo = u_norm_periodic(xi_t - eps, w0_lh)
    u_hi = u_norm_periodic(xi_t + eps, w0_lh)
    local_xit_exp = math.log(u_hi / u_lo) / math.log((xi_t + eps) / (xi_t - eps))

    # Locate exact maximum of u_norm(xi_t)
    res = minimize_scalar(
        lambda x: -u_norm_periodic(x, w0_lh),
        bounds=(5.0, 500.0),
        method="bounded",
        options={"xatol": opt_tol},
    )
    xi_t_max = float(res.x)
    beta_U_crit = math.sqrt(4.0 / u_n) if u_n > 0.0 else float("nan")

    return {
        "u": u_val,
        "u_norm": u_n,
        "xi_t_max": xi_t_max,
        "stiffness_exponent": stiffness_exp,
        "local_xit_exponent": local_xit_exp,
        "beta_U_crit": beta_U_crit,
    }
