"""Eq. (2) evaluated as a function of L^2 / delta for Fig. 4."""

import math


def run(inputs, settings):
    kappa_f = inputs["kappa_f"]
    N = inputs["N"]
    alpha = (inputs["k_x"] / inputs["E_f_A_f"]) * inputs["L2_over_delta"]
    chi2 = inputs["chi2"]
    c = inputs["c"]
    kappa_B = kappa_f * N * (1.0 + chi2 * (N - 1.0) / (1.0 + c * (N + math.sqrt(N)) / alpha))
    return {"kappa_B": kappa_B}


def symbolic():
    import sympy as sp

    kappa_f, N, k_x, E_f_A_f, L2_over_delta, chi2, c = sp.symbols(
        "kappa_f N k_x E_f_A_f L2_over_delta chi2 c", positive=True
    )
    alpha = (k_x / E_f_A_f) * L2_over_delta
    return {"kappa_B": kappa_f * N * (1 + chi2 * (N - 1) / (1 + c * (N + sp.sqrt(N)) / alpha))}
