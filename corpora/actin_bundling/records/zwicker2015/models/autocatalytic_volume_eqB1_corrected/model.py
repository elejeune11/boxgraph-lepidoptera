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

    V0 = (V_s / psi) * (phibar - k_b / k)
    delta_V = (Q * k_b - 4.0 * math.pi * D_B * k * gb) / ((k * phibar - k_b) * k * psi)
    return {"V_bar": V0 + delta_V, "delta_V": delta_V}


def symbolic():
    import sympy as sp

    V_s, psi, phibar, k_b, k, Q, D_B, gamma_beta = sp.symbols(
        "V_s psi phibar k_b k Q D_B gamma_beta", positive=True
    )
    V0 = (V_s / psi) * (phibar - k_b / k)
    delta_V = (Q * k_b - 4 * sp.pi * D_B * k * gamma_beta) / ((k * phibar - k_b) * k * psi)
    return {"V_bar": V0 + delta_V, "delta_V": delta_V}
