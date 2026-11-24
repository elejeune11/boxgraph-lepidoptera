import math


def run(inputs, settings):
    phi_0_A = inputs["k_b"] / inputs["k"]
    return {"Q_c": 4.0 * math.pi * inputs["D_B"] * inputs["gamma_beta"] / phi_0_A}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, k_b, k = sp.symbols("D_B gamma_beta k_b k", positive=True)
    return {"Q_c": 4 * sp.pi * D_B * gamma_beta * k / k_b}
