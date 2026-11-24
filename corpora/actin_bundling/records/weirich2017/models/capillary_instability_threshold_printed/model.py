import math


def run(inputs, settings):
    enh = math.sqrt(1.0 + inputs["omega"])
    lam_c_um = 2.0 * math.pi * (inputs["R_0"] * 1e-3) * enh
    return {"lambda_c": lam_c_um, "enhancement": enh}


def symbolic():
    import sympy as sp

    R_0, omega = sp.symbols("R_0 omega", positive=True)
    enh = sp.sqrt(1 + omega)
    return {
        "lambda_c": 2 * sp.pi * sp.Rational(1, 1000) * R_0 * enh,
        "enhancement": enh,
    }
