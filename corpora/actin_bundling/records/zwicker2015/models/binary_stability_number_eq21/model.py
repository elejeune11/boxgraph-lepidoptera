import math


def run(inputs, settings):
    V_s = inputs["V_s"]
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    N = inputs["N"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    n_c = (1.0 / (2.0 * math.pi * D_B * gb)) * (k_f * k_b / (k_f + k_b))
    N_max = V_s * n_c
    chi = N_max / N
    return {"chi": chi, "N_max": N_max, "n_c": n_c}


def symbolic():
    import sympy as sp

    V_s, D_B, gamma_beta, N, k_f, k_b = sp.symbols("V_s D_B gamma_beta N k_f k_b", positive=True)
    n_c = (1 / (2 * sp.pi * D_B * gamma_beta)) * (k_f * k_b / (k_f + k_b))
    N_max = V_s * n_c
    chi = N_max / N
    return {"chi": chi, "N_max": N_max, "n_c": n_c}
