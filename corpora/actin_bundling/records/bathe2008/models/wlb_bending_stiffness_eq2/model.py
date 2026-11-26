"""Closed-form Wormlike Bundle (WLB) bending stiffness: Eq. (2) and Eq. (A6)."""

import math


def run(inputs, settings):
    kappa_f = inputs["kappa_f"]
    N = inputs["N"]
    alpha = inputs["alpha"]
    chi2 = inputs["chi2"]
    c = inputs["c"]
    kappa_B_star = N * (1.0 + chi2 * (N - 1.0) / (1.0 + c * (N + math.sqrt(N)) / alpha))
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}


def symbolic():
    import sympy as sp

    kappa_f, N, alpha, chi2, c = sp.symbols("kappa_f N alpha chi2 c", positive=True)
    kappa_B_star = N * (1 + chi2 * (N - 1) / (1 + c * (N + sp.sqrt(N)) / alpha))
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}
