import math


def _bracket(beta, L_B):
    u = 0.5 * beta * L_B
    if u < 1e-3:
        series = (
            (2.0 / 15.0) * (u**5)
            - (17.0 / 315.0) * (u**7)
            + (62.0 / 2835.0) * (u**9)
        )
        return series / (beta**5)
    beta2 = beta * beta
    return (L_B**3) / (24.0 * beta2) - L_B / (2.0 * beta2 * beta2) + math.tanh(u) / (beta**5)


def run(inputs, settings):
    n3d = inputs["n"]
    L_B = inputs["L_B"]
    k_par = inputs["k_par"]
    delta = inputs["delta"]
    k_a = inputs["k_a"]
    kappa_a = inputs["kappa_a"]
    d = inputs["d"]
    t = inputs["t"]

    n2d = math.sqrt(n3d)
    if n2d <= 1.0 or k_par <= 0.0:
        return {"kappa_B": n3d * kappa_a}

    S2 = n2d * (n2d * n2d - 1.0) / 24.0
    K = n2d * (n2d + 1.0) / 12.0
    geom = d * (d + t)

    beta2 = (k_par / (delta * K * k_a)) + (2.0 * k_par * geom * S2) / (delta * K * n2d * kappa_a)
    beta = math.sqrt(beta2)
    br = _bracket(beta, L_B)
    w_max = (L_B**3) / (48.0 * n2d * kappa_a) - (
        k_par * geom * S2 / (delta * K * (n2d**2) * (kappa_a**2))
    ) * br
    kappa_B_2d = (L_B**3) / (48.0 * w_max)
    return {"kappa_B": n2d * kappa_B_2d}


def symbolic():
    import sympy as sp

    n, L_B, k_par, delta, k_a, kappa_a, d, t = sp.symbols(
        "n L_B k_par delta k_a kappa_a d t", positive=True
    )
    n2d = sp.sqrt(n)
    S2 = n2d * (n2d**2 - 1) / 24
    K = n2d * (n2d + 1) / 12
    geom = d * (d + t)
    beta2 = k_par / (delta * K * k_a) + 2 * k_par * geom * S2 / (delta * K * n2d * kappa_a)
    beta = sp.sqrt(beta2)
    br = L_B**3 / (24 * beta2) - L_B / (2 * beta2**2) + sp.tanh(beta * L_B / 2) / (beta**5)
    w_max = L_B**3 / (48 * n2d * kappa_a) - (
        k_par * geom * S2 / (delta * K * n2d**2 * kappa_a**2)
    ) * br
    kappa_B_2d = L_B**3 / (48 * w_max)
    return {"kappa_B": n2d * kappa_B_2d}
