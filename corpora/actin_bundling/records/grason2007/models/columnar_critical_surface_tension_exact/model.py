import math
from scipy.optimize import brentq


def run(inputs, settings):
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    delta_mu = inputs["delta_mu"]
    rho_0 = inputs["rho_0"]
    L = inputs["L"]
    rtol = float(settings.get("rtol", 1e-8))

    if delta_mu <= 0.0:
        u_star = 1.0
        h_val = 2.0 ** (-1.0 / 3.0)
    else:
        d_tilde = math.pi * (rho_0 / L) * delta_mu * (lam3 ** (4.0 / 3.0)) / f_0

        def hp(u):
            return (1.0 - u**2) / (3.0 * (u ** (2.0 / 3.0)) * ((1.0 + u**2) ** (4.0 / 3.0))) - d_tilde

        # Asymptotic root for large d_tilde is u ~ (3 * d_tilde)^(-3/2)
        u_lo = min(0.5, 0.25 * (3.0 * d_tilde) ** (-1.5))
        while hp(u_lo) <= 0.0:
            u_lo *= 0.25
        u_star = brentq(hp, u_lo, 1.0, rtol=rtol, xtol=1e-300)
        h_val = (u_star ** (1.0 / 3.0)) / ((1.0 + u_star**2) ** (1.0 / 3.0)) - d_tilde * u_star

    sigma_c = f_0 * h_val / (2.0 * math.pi * (lam3 ** (1.0 / 3.0)))
    return {"Sigma_c": sigma_c, "R_star": u_star * lam3}
