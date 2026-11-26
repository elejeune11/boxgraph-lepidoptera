import math
from scipy.special import k0, k1


def run(inputs, settings):
    D = inputs["D"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    l_B = inputs["l_B"]
    lam = inputs["lambda_e"]
    kappa = 1.0 / kappa_inv

    x = kappa * D / 2.0
    denom = x * float(k1(x))
    pref = (lam ** 2) * l_B

    Gamma_par = 2.0 * D * pref * float(k0(kappa * (D + delta))) / (denom ** 2)
    Gamma_perp = 2.0 * math.pi * kappa_inv * pref * math.exp(-kappa * (D + delta)) / (denom ** 2)
    return {"Gamma_par": float(Gamma_par), "Gamma_perp": float(Gamma_perp)}


def symbolic():
    import sympy as sp

    D, delta, kappa_inv, l_B, lambda_e = sp.symbols("D delta kappa_inv l_B lambda_e", positive=True)
    kappa = 1 / kappa_inv
    x = kappa * D / 2
    denom = x * sp.besselk(1, x)
    pref = lambda_e ** 2 * l_B
    return {
        "Gamma_par": 2 * D * pref * sp.besselk(0, kappa * (D + delta)) / denom ** 2,
        "Gamma_perp": 2 * sp.pi * kappa_inv * pref * sp.exp(-kappa * (D + delta)) / denom ** 2,
    }
