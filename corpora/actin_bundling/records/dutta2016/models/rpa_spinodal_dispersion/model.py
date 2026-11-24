import math


def _debye_b(k):
    k = abs(k)
    if k < 0.1:
        k2 = k * k
        return (
            1.0
            - k2 / 6.0
            + (k2**2) / 48.0
            - (k2**3) / 480.0
            + (k2**4) / 5760.0
            - (k2**5) / 80640.0
        )
    return (-8.0 + 8.0 * math.exp(-0.5 * k * k) + 4.0 * k * k) / (k**4)


def run(inputs, settings):
    k = float(inputs["k"])
    rho_star = float(inputs["rho_star"])
    GR = float(inputs["Gamma_R"])
    lR = float(inputs["lambda_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    bv = (
        2.0 * math.pi * GR / (k * k + 1.0 / (lR * lR))
        - math.pi * GA * (lA * lA) * math.exp(-0.25 * (k * k) * (lA * lA))
    )
    bk = _debye_b(k)
    g_static = -rho_star * bv * bk - 1.0
    lam_B = -(k * k) * (1.0 / (rho_star * bk) + bv)
    return {"g_static": g_static, "lam_B": lam_B}


def symbolic():
    import sympy as sp

    k, rho_star, Gamma_R, lambda_R, Gamma_A, lambda_A = sp.symbols(
        "k rho_star Gamma_R lambda_R Gamma_A lambda_A", positive=True
    )
    bv = (
        2 * sp.pi * Gamma_R / (k**2 + 1 / lambda_R**2)
        - sp.pi * Gamma_A * lambda_A**2 * sp.exp(-(k**2) * lambda_A**2 / 4)
    )
    bk = (-8 + 8 * sp.exp(-(k**2) / 2) + 4 * k**2) / k**4
    g_static = -rho_star * bv * bk - 1
    lam_B = -(k**2) * (1 / (rho_star * bk) + bv)
    return {"g_static": g_static, "lam_B": lam_B}
