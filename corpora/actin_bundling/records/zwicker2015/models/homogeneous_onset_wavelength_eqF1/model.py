import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    w = inputs["w"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]

    ksum = k_f + k_b
    q_c = (4.0 * ksum / (D_B * w * w)) ** 0.25
    inv_qc = 1.0 / q_c
    lambda_c = 2.0 * math.pi * inv_qc
    return {"q_c": q_c, "inv_qc": inv_qc, "lambda_c": lambda_c}


def symbolic():
    import sympy as sp

    D_B, w, k_f, k_b = sp.symbols("D_B w k_f k_b", positive=True)
    q_c = (4 * (k_f + k_b) / (D_B * w**2)) ** sp.Rational(1, 4)
    inv_qc = (D_B * w**2 / (4 * (k_f + k_b))) ** sp.Rational(1, 4)
    lambda_c = 2 * sp.pi * inv_qc
    return {"q_c": q_c, "inv_qc": inv_qc, "lambda_c": lambda_c}
