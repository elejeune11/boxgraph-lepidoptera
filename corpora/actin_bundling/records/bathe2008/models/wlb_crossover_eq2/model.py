"""Factor-of-2 regime crossover values alpha_dec and alpha_fc solved from Eq. (2)."""

import math


def run(inputs, settings):
    N = inputs["N"]
    chi2 = inputs["chi2"]
    c = inputs["c"]
    A = chi2 * (N - 1.0)
    # Decoupled-to-intermediate: kappa_B* = 2*N => 1 + A/(1 + c*(N+sqrt(N))/alpha_dec) = 2
    alpha_dec = c * (N + math.sqrt(N)) / (A - 1.0)
    # Fully coupled-to-intermediate: kappa_B* = (1/2)*N*(1 + A) => 1 + A/(1 + c*(N+sqrt(N))/alpha_fc) = (1+A)/2
    alpha_fc = c * (N + math.sqrt(N)) * (A - 1.0) / (A + 1.0)
    return {"alpha_dec": alpha_dec, "alpha_fc": alpha_fc}


def symbolic():
    import sympy as sp

    N, chi2, c = sp.symbols("N chi2 c", positive=True)
    A = chi2 * (N - 1)
    return {
        "alpha_dec": c * (N + sp.sqrt(N)) / (A - 1),
        "alpha_fc": c * (N + sp.sqrt(N)) * (A - 1) / (A + 1),
    }
