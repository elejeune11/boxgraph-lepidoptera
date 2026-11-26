"""Microscopic-mapped Timoshenko bundle bending stiffness: kappa_B = N^2 * kappa_f / (1 + N/alpha)."""


def run(inputs, settings):
    kappa_f = inputs["kappa_f"]
    N = inputs["N"]
    alpha = inputs["alpha"]
    kappa_B_star = (N ** 2) / (1.0 + N / alpha)
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}


def symbolic():
    import sympy as sp

    kappa_f, N, alpha = sp.symbols("kappa_f N alpha", positive=True)
    kappa_B_star = N**2 / (1 + N / alpha)
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}
