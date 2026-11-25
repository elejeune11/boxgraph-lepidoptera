"""Linear stability analysis of the exact mirza2025_code continuum PDEs.

Ports the linearization of AuxCortexNematic2D.cpp lines 420, 514-565, 1070-1140
around the isotropic uniform steady state (h = heqb, q = 0, v = 0), including:
- L_eff = frank / 2 (from DqDq = 0.5*(|grad q1|^2 + |grad q2|^2))
- c0_eff = 2*sus20 + lambda_rot*heqb + 4*k_lim/S2_max (including the hidden log barrier at q=0)
- kappa = (lambda_aniso - crossLinker) / lambda_iso
Computes marginal stability (lambda_crit, nu_crit, ell_pattern), checks config's
lambda_iso / lambda_crit (= 1.3 to 16 digits at k_lim=0) and r_out / ell_pattern (= 8.043),
and computes the fastest-growing supercritical mode ratio nu_crit / nu_max(lambda_iso).
"""

import math
import numpy as np
from scipy.optimize import minimize_scalar


def run(inputs, settings):
    kd = float(inputs["visc"])
    beta = float(inputs["cvisc"])
    er = float(inputs["rvisc"])
    L_eff = 0.5 * float(inputs["frank"])
    sus20 = float(inputs["sus20"])
    lam_rot = float(inputs["lambda_rot"])
    lam_iso = float(inputs["lambda_iso"])
    lam_aniso = float(inputs["lambda_aniso"])
    xl = float(inputs["crossLinker"])
    heqb = float(inputs["heqb"])
    k_lim = float(inputs["k_lim"])
    S2_max = float(inputs["S2_max"])
    r_out = float(inputs["r_out"])

    c0_eff = 2.0 * sus20 + lam_rot * heqb + 4.0 * k_lim / S2_max
    kap = (lam_aniso - xl) / lam_iso

    def lam_marginal(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + kd) * (2.0 * c0_eff + 2.0 * L_eff * x)
        den = (2.0 * c0_eff + 2.0 * L_eff * x) + kap * beta * (x + kd)
        return num / den if den > 0.0 else 1e30

    n_grid = int(round(float(settings.get("n_grid", 4000))))
    x_min = float(settings.get("x_min", 1e-4))
    x_max = float(settings.get("x_max", 1e3))

    xs = np.geomspace(x_min, x_max, n_grid)
    vals = np.array([lam_marginal(float(x)) for x in xs])
    i0 = int(np.argmin(vals))
    res = minimize_scalar(
        lam_marginal,
        bounds=(float(xs[max(0, i0 - 3)]), float(xs[min(n_grid - 1, i0 + 3)])),
        method="bounded",
        options={"xatol": 1e-13},
    )
    x_c = float(res.x)
    lam_c = float(res.fun)
    nu_c = math.sqrt(x_c)
    ell_pat = 2.0 * math.pi / nu_c

    # Supercritical fastest-growing linear mode at lambda = lam_iso:
    # Quadratic dispersion relation A(x)*alpha^2 + B(x)*alpha + C(x) = 0
    # matching mirza2025 Eq. (C7) with L_eff and c0_eff at q0 = 0
    def max_growth(x):
        if x <= 0.0:
            return -1e30
        A = (2.0 / x + 8.0) * er - beta * beta
        B = (
            (2.0 / x + 8.0) * (c0_eff + L_eff * x)
            + (x + kd) * (2.0 * er / x + 8.0 * er - beta * beta)
            - 2.0 * er * lam_iso
        )
        C = (
            (x + kd) * (2.0 / x + 8.0) * (c0_eff + L_eff * x)
            - 2.0 * lam_iso * (c0_eff + L_eff * x)
            - kap * beta * lam_iso * (x + kd)
        )
        disc = B * B - 4.0 * A * C
        if disc < 0.0 or A <= 0.0:
            return -1e30
        return (-B + math.sqrt(disc)) / (2.0 * A)

    gvals = np.array([max_growth(float(x)) for x in xs])
    ig = int(np.argmax(gvals))
    res_g = minimize_scalar(
        lambda x: -max_growth(x),
        bounds=(float(xs[max(0, ig - 3)]), float(xs[min(n_grid - 1, ig + 3)])),
        method="bounded",
        options={"xatol": 1e-13},
    )
    nu_max = math.sqrt(float(res_g.x))

    return {
        "lambda_crit": lam_c,
        "nu_crit": nu_c,
        "ell_pattern": ell_pat,
        "activity_ratio": lam_iso / lam_c,
        "domain_wavelengths": r_out / ell_pat,
        "nu_max_supercrit": nu_max,
        "wavelength_ratio_lin": nu_c / nu_max,
    }
