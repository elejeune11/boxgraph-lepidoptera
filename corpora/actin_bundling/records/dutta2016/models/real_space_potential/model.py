import math
from scipy.special import k0


def run(inputs, settings):
    r = inputs["r"]
    vr = inputs["Gamma_R"] * float(k0(r / inputs["lambda_R"]))
    va = -inputs["Gamma_A"] * math.exp(-(r ** 2) / (inputs["lambda_A"] ** 2))
    return {"betaV": vr + va}


def symbolic():
    import sympy as sp
    r, Gamma_R, lambda_R, Gamma_A, lambda_A = sp.symbols(
        "r Gamma_R lambda_R Gamma_A lambda_A", positive=True
    )
    return {
        "betaV": Gamma_R * sp.besselk(0, r / lambda_R)
        - Gamma_A * sp.exp(-(r ** 2) / (lambda_A ** 2))
    }
