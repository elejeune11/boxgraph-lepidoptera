import math


def run(inputs, settings):
    kappa = inputs["kappa"]
    H = inputs["H"]
    n = inputs["n"]
    H0 = inputs["H0"]
    sigma = inputs["sigma"]
    kBT = inputs["kBT"]
    ns = inputs["ns"]
    J = inputs["J"]
    grad_n = inputs["grad_n"]

    f_bend = 0.5 * kappa * (H - n * H0) ** 2
    f_ent = kBT * ns * n * (math.log(n) - 1.0)
    f_int = -J * (n ** 2) + (J / ns) * (grad_n ** 2)
    return {"f": f_bend + sigma + f_ent + f_int}


def symbolic():
    import sympy as sp

    kappa, H, n, H0, sigma, kBT, ns, J, grad_n = sp.symbols(
        "kappa H n H0 sigma kBT ns J grad_n", real=True
    )
    f = (
        sp.Rational(1, 2) * kappa * (H - n * H0) ** 2
        + sigma
        + kBT * ns * n * (sp.log(n) - 1)
        - J * n ** 2
        + (J / ns) * grad_n ** 2
    )
    return {"f": f}
