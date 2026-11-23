import math


def _debye_b(k):
    k = abs(k)
    if k < 1e-2:
        k2 = k * k
        return 1.0 - k2 / 6.0 + (k2 * k2) / 48.0 - (k2 * k2 * k2) / 480.0
    return (-8.0 + 8.0 * math.exp(-0.5 * k * k) + 4.0 * k * k) / (k ** 4)


def run(inputs, settings):
    k = float(inputs["k"])
    q = 2.0 * math.pi * k
    rho_star = float(inputs["rho_star"])
    GR = float(inputs["Gamma_R"])
    lR = float(inputs["lambda_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    bv = 2.0 * math.pi * GR / (q * q + 1.0 / (lR * lR)) - math.pi * GA * (lA * lA) * math.exp(-0.25 * (q * q) * (lA * lA))
    return {"F": rho_star * abs(bv) * _debye_b(k)}


def symbolic():
    import sympy as sp
    k, rho_star, Gamma_R, lambda_R, Gamma_A, lambda_A = sp.symbols(
        "k rho_star Gamma_R lambda_R Gamma_A lambda_A", positive=True
    )
    q = 2 * sp.pi * k
    bv = 2 * sp.pi * Gamma_R / (q**2 + 1 / lambda_R**2) - sp.pi * Gamma_A * lambda_A**2 * sp.exp(-(q**2) * lambda_A**2 / 4)
    bk = (-8 + 8 * sp.exp(-k**2 / 2) + 4 * k**2) / k**4
    return {"F": rho_star * sp.Abs(bv) * bk}
