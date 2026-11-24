import math


def run(inputs, settings):
    k_max = 1.0 / math.sqrt(2.0 * (1.0 + 2.0 * inputs["omega"]))
    ratio = 2.0 * math.pi / k_max
    lam_max_um = ratio * (inputs["R_0"] * 1e-3)
    return {
        "k_max": k_max,
        "lambda_max": lam_max_um,
        "lambda_to_R0": ratio,
    }


def symbolic():
    import sympy as sp

    R_0, omega = sp.symbols("R_0 omega", positive=True)
    k_max = 1 / sp.sqrt(2 * (1 + 2 * omega))
    ratio = 2 * sp.pi * sp.sqrt(2 * (1 + 2 * omega))
    return {
        "k_max": k_max,
        "lambda_max": ratio * sp.Rational(1, 1000) * R_0,
        "lambda_to_R0": ratio,
    }
