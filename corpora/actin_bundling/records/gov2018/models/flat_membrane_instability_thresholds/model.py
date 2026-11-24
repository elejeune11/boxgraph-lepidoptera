import math


def run(inputs, settings):
    kappa = inputs["kappa"]
    sigma = inputs["sigma"]
    kBT = inputs["kBT"]
    H0_mag = abs(inputs["H0"])
    ns = inputs["ns"]
    n0 = inputs["n0"]
    A_mag = abs(inputs["A"])

    A_c = sigma * H0_mag + sigma * kBT * ns / (kappa * H0_mag * n0)
    sigma_c = (kappa * H0_mag * n0 * A_mag) / (kBT * ns + kappa * (H0_mag ** 2) * n0)
    num = kappa * H0_mag * n0 * (A_mag - sigma * H0_mag) - sigma * kBT * ns
    if num > 0.0:
        q_c = math.sqrt(num / (kappa * kBT * ns))
        lambda_c = 2.0 * math.pi / q_c
    else:
        q_c = None
        lambda_c = None
    return {"A_c": A_c, "sigma_c": sigma_c, "q_c": q_c, "lambda_c": lambda_c}


def symbolic():
    import sympy as sp

    kappa, sigma, kBT, H0, ns, n0, A = sp.symbols(
        "kappa sigma kBT H0 ns n0 A", positive=True
    )
    A_c = sigma * H0 + sigma * kBT * ns / (kappa * H0 * n0)
    sigma_c = (kappa * H0 * n0 * A) / (kBT * ns + kappa * H0 ** 2 * n0)
    q_c = sp.sqrt((kappa * H0 * n0 * (A - sigma * H0) - sigma * kBT * ns) / (kappa * kBT * ns))
    return {
        "A_c": A_c,
        "sigma_c": sigma_c,
        "q_c": q_c,
        "lambda_c": 2 * sp.pi / q_c,
    }
