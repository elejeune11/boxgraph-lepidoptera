import math
from scipy.optimize import brentq, minimize_scalar


def _omega0_eq9(R, mu_par_bar, mu_perp_bar, K3_bar, gamma, rtol):
    a = 4.0 * (mu_perp_bar * (R**4) + K3_bar * (R**2))
    b = 2.0 * mu_par_bar * (R**2)
    hi = min(gamma / b, (gamma / a) ** (1.0 / 3.0)) if b > 0.0 else (gamma / a) ** (1.0 / 3.0)
    return brentq(
        lambda om: a * (om**3) + b * om - gamma,
        0.0,
        hi * (1.0 + 1e-12),
        rtol=rtol,
        xtol=1e-300,
    )


def run(inputs, settings):
    mu_par_bar = inputs["mu_par_bar"]
    mu_perp_bar = inputs["mu_perp_bar"]
    K3_bar = inputs["K3_bar"]
    gamma = inputs["gamma"]
    Sigma = inputs["Sigma"]
    rtol = float(settings.get("rtol", 1e-8))
    inner_rtol = min(rtol * 1e-2, 1e-13)

    # By the envelope theorem on Eq. (9): d/dR [pi*R^2*E_elast_min(R) + 2*pi*R*Sigma] = 0
    # <=> lhs(R) = Omega_0(R)^4 * R^3 * (mu_perp_bar*R^2 + 2*K3_bar) = Sigma.
    def lhs(R):
        om = _omega0_eq9(R, mu_par_bar, mu_perp_bar, K3_bar, gamma, inner_rtol)
        return (om**4) * (R**3) * (mu_perp_bar * (R**2) + 2.0 * K3_bar)

    lam3 = math.sqrt(K3_bar / mu_perp_bar)
    r_char = gamma / mu_par_bar
    log_lo = math.log(min(lam3, r_char)) - 18.0
    log_hi = math.log(max(lam3, r_char)) + 18.0
    res = minimize_scalar(
        lambda u: -lhs(math.exp(u)),
        bounds=(log_lo, log_hi),
        method="bounded",
        options={"xatol": min(rtol, 1e-8)},
    )
    r_peak = math.exp(res.x)
    peak_val = -res.fun
    if Sigma > peak_val * (1.0 + 1e-12):
        return {"R_0": None}
    if abs(Sigma - peak_val) <= 1e-12 * peak_val:
        return {"R_0": r_peak}

    r_hi = r_peak * 2.0
    while lhs(r_hi) > Sigma:
        r_hi *= 2.0
    r_0 = brentq(lambda R: lhs(R) - Sigma, r_peak, r_hi, rtol=rtol, xtol=1e-300)
    return {"R_0": r_0}
