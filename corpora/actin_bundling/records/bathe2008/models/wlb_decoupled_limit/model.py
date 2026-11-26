"""Decoupled bundle limit (alpha << 1): kappa_B = N * kappa_f."""


def run(inputs, settings):
    N = inputs["N"]
    return {"kappa_B": N * inputs["kappa_f"], "kappa_B_star": N}


def symbolic():
    import sympy as sp

    kappa_f, N = sp.symbols("kappa_f N", positive=True)
    return {"kappa_B": N * kappa_f, "kappa_B_star": N}
