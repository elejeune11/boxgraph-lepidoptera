"""Fully coupled bundle limit (alpha >> N): kappa_B = kappa_f * N * (1 + chi^2 * (N - 1))."""


def run(inputs, settings):
    N = inputs["N"]
    chi2 = inputs["chi2"]
    kb_star = N * (1.0 + chi2 * (N - 1.0))
    return {"kappa_B": inputs["kappa_f"] * kb_star, "kappa_B_star": kb_star}


def symbolic():
    import sympy as sp

    kappa_f, N, chi2 = sp.symbols("kappa_f N chi2", positive=True)
    kb_star = N * (1 + chi2 * (N - 1))
    return {"kappa_B": kappa_f * kb_star, "kappa_B_star": kb_star}
