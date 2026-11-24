import math


def run(inputs, settings):
    V_s = inputs["V_s"]
    psi = inputs["psi"]
    phibar = inputs["phibar"]
    k_b = inputs["k_b"]
    k = inputs["k"]
    Q = inputs["Q"]
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]

    A = k * psi**2
    B = -V_s * psi * (k * phibar - k_b - Q / V_s)
    C = -V_s * (Q * phibar - 4.0 * math.pi * D_B * gb)
    disc = B**2 - 4.0 * A * C
    if disc < 0.0:
        return {"V_bar": None, "delta_V": None}
    V_bar = (-B + math.sqrt(disc)) / (2.0 * A)
    V0 = (V_s / psi) * (phibar - k_b / k)
    return {"V_bar": V_bar, "delta_V": V_bar - V0}


def symbolic():
    import sympy as sp

    V_s, psi, phibar, k_b, k, Q, D_B, gamma_beta = sp.symbols(
        "V_s psi phibar k_b k Q D_B gamma_beta", positive=True
    )
    A = k * psi**2
    B = -V_s * psi * (k * phibar - k_b - Q / V_s)
    C = -V_s * (Q * phibar - 4 * sp.pi * D_B * gamma_beta)
    V_bar = (-B + sp.sqrt(B**2 - 4 * A * C)) / (2 * A)
    V0 = (V_s / psi) * (phibar - k_b / k)
    return {"V_bar": V_bar, "delta_V": V_bar - V0}
