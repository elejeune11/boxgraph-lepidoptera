import cmath
import math
from scipy.optimize import minimize_scalar

# Printed six cross-link sites per 28-monomer unit cell (p. 3 col 1 & Fig. 2c)
SITES_28_13 = ((0, 0), (4, 1), (9, -1), (14, 3), (19, 1), (24, -1))


def u_norm_periodic(xi_t: float, w0: float) -> float:
    """Exact infinite-chain periodic geometric sum for u / (beta*U)^2 from Eq. (3)."""
    z_p = cmath.exp(-28.0 / xi_t + 1j * 28.0 * w0)
    total = 0.0
    for l1, m1 in SITES_28_13:
        for l2, m2 in SITES_28_13:
            dphi = 2.0 * math.pi * (m1 - m2) / 6.0
            dl0 = l1 - l2
            # Positive shifts: l1 - (l2 + 28*k) > 0
            d_pos = dl0 % 28
            if d_pos == 0:
                d_pos = 28
            s_pos = (cmath.exp(-d_pos / xi_t + 1j * (w0 * d_pos - dphi)) / (1.0 - z_p)).real
            # Negative shifts: l1 - (l2 + 28*k) < 0
            d_neg = (-dl0) % 28
            if d_neg == 0:
                d_neg = 28
            s_neg = (cmath.exp(-d_neg / xi_t + 1j * (w0 * d_neg + dphi)) / (1.0 - z_p)).real
            total += 0.5 * (s_pos + s_neg)
    return total / float(len(SITES_28_13))


def run(inputs, settings):
    beta_U = float(inputs["beta_U"])
    xi_t = float(inputs["xi_t"])
    w0 = float(inputs["omega_0"])
    opt_tol = float(settings.get("opt_tol", 1.0e-6))

    u_n = u_norm_periodic(xi_t, w0)
    u_val = (beta_U ** 2) * u_n

    # Log-slope with respect to beta_U (exact 2 for nonzero u_n)
    b1, b2 = beta_U * 0.99, beta_U * 1.01
    u1, u2 = (b1 ** 2) * u_n, (b2 ** 2) * u_n
    stiffness_exp = math.log(abs(u2 / u1)) / math.log(b2 / b1) if abs(u_n) > 0 else float("nan")

    # Maximize u_norm(xi_t) over [1, 1000]
    res = minimize_scalar(
        lambda x: -u_norm_periodic(x, w0),
        bounds=(1.0, 1000.0),
        method="bounded",
        options={"xatol": opt_tol},
    )
    xi_t_max = float(res.x)
    u_max = float(-res.fun)

    return {
        "u": u_val,
        "u_norm": u_n,
        "xi_t_max": xi_t_max,
        "u_max": u_max,
        "stiffness_exponent": stiffness_exp,
    }
