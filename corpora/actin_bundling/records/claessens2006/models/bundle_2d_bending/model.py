import math


def _bracket(beta, L_B):
    """Evaluate L_B^3/(24*beta^2) - L_B/(2*beta^4) + tanh(beta*L_B/2)/beta^5
    without catastrophic cancellation for small beta*L_B."""
    u = 0.5 * beta * L_B
    if u < 1e-3:
        # u^3/3 - u + tanh(u) = 2*u^5/15 - 17*u^7/315 + 62*u^9/2835 - ...
        series = (
            (2.0 / 15.0) * (u**5)
            - (17.0 / 315.0) * (u**7)
            + (62.0 / 2835.0) * (u**9)
        )
        return series / (beta**5)
    beta2 = beta * beta
    return (L_B**3) / (24.0 * beta2) - L_B / (2.0 * beta2 * beta2) + math.tanh(u) / (beta**5)


def run(inputs, settings):
    n = inputs["n"]
    L_B = inputs["L_B"]
    k_par = inputs["k_par"]
    delta = inputs["delta"]
    k_a = inputs["k_a"]
    kappa_a = inputs["kappa_a"]
    d = inputs["d"]
    t = inputs["t"]

    N = 0.5 * (n - 1.0)
    # Exact closed forms for integer N >= 1 (and smooth continuation for n > 1):
    # sum_{i=1}^N i^2 = N*(N+1)*(2N+1)/6 = n*(n^2 - 1)/24
    # K = (1/N)*sum_{k=1}^N sum_{i=k}^N i = (1/N)*sum_{i=1}^N i^2 = n*(n+1)/12
    S2 = n * (n * n - 1.0) / 24.0
    K = n * (n + 1.0) / 12.0
    geom = d * (d + t)

    if n <= 1.0 or k_par <= 0.0:
        w_max = (L_B**3) / (48.0 * n * kappa_a)
        beta = math.sqrt(max(0.0, k_par / (delta * K * k_a))) if k_par > 0.0 else 0.0
        return {"beta": beta, "w_B_max": w_max, "kappa_B_2d": n * kappa_a}

    beta2 = (k_par / (delta * K * k_a)) + (2.0 * k_par * geom * S2) / (delta * K * n * kappa_a)
    beta = math.sqrt(beta2)
    br = _bracket(beta, L_B)
    w_max = (L_B**3) / (48.0 * n * kappa_a) - (
        k_par * geom * S2 / (delta * K * (n**2) * (kappa_a**2))
    ) * br
    kappa_B_2d = (L_B**3) / (48.0 * w_max)
    return {"beta": beta, "w_B_max": w_max, "kappa_B_2d": kappa_B_2d}


def symbolic():
    import sympy as sp

    n, L_B, k_par, delta, k_a, kappa_a, d, t = sp.symbols(
        "n L_B k_par delta k_a kappa_a d t", positive=True
    )
    S2 = n * (n**2 - 1) / 24
    K = n * (n + 1) / 12
    geom = d * (d + t)
    beta2 = k_par / (delta * K * k_a) + 2 * k_par * geom * S2 / (delta * K * n * kappa_a)
    beta = sp.sqrt(beta2)
    br = L_B**3 / (24 * beta2) - L_B / (2 * beta2**2) + sp.tanh(beta * L_B / 2) / (beta**5)
    w_max = L_B**3 / (48 * n * kappa_a) - (
        k_par * geom * S2 / (delta * K * n**2 * kappa_a**2)
    ) * br
    kappa_B_2d = L_B**3 / (48 * w_max)
    return {"beta": beta, "w_B_max": w_max, "kappa_B_2d": kappa_B_2d}
