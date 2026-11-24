import math
from scipy.special import lambertw


def run(inputs, settings):
    kappa = inputs["kappa"]
    kBT = inputs["kBT"]
    H = inputs["H"]
    H0 = inputs["H0"]
    ns = inputs["ns"]
    n0 = inputs["n0"]

    alpha = (kappa * (H0 ** 2)) / (kBT * ns)
    chi = (kappa * H * H0) / (kBT * ns)
    z = alpha * n0 * math.exp(alpha * n0 + chi)
    n = float(lambertw(z).real) / alpha
    return {"enrichment": n / n0}


def symbolic():
    import sympy as sp

    kappa, kBT, H, H0, ns, n0 = sp.symbols("kappa kBT H H0 ns n0", positive=True)
    alpha = (kappa * H0 ** 2) / (kBT * ns)
    chi = (kappa * H * H0) / (kBT * ns)
    return {"enrichment": sp.LambertW(alpha * n0 * sp.exp(alpha * n0 + chi)) / (alpha * n0)}
