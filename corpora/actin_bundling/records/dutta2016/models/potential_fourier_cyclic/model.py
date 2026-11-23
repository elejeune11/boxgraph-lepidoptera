import math


def run(inputs, settings):
    q = 2.0 * math.pi * float(inputs["k"])
    GR = float(inputs["Gamma_R"])
    lR = float(inputs["lambda_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    vr_k = 2.0 * math.pi * GR / (q * q + 1.0 / (lR * lR))
    va_k = -math.pi * GA * (lA * lA) * math.exp(-0.25 * (q * q) * (lA * lA))
    return {"betaV": vr_k + va_k}


def symbolic():
    import sympy as sp
    k, Gamma_R, lambda_R, Gamma_A, lambda_A = sp.symbols(
        "k Gamma_R lambda_R Gamma_A lambda_A", positive=True
    )
    q = 2 * sp.pi * k
    vr_k = 2 * sp.pi * Gamma_R / (q**2 + 1 / lambda_R**2)
    va_k = -sp.pi * Gamma_A * lambda_A**2 * sp.exp(-(q**2) * lambda_A**2 / 4)
    return {"betaV": vr_k + va_k}
