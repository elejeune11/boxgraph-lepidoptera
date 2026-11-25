import math

NA_UM3_PER_UMOL = 602.214076  # molecules per (uM * um^3)


def run(inputs, settings):
    R = float(inputs["R_drop"])
    V_drop = (4.0 / 3.0) * math.pi * (R ** 3)
    C = float(inputs["N"]) / (NA_UM3_PER_UMOL * V_drop)
    return {"C": C, "V_drop": V_drop}


def symbolic():
    import sympy as sp

    N, R_drop = sp.symbols("N R_drop", positive=True)
    V_drop = sp.Rational(4, 3) * sp.pi * R_drop**3
    C = N / (sp.Float(NA_UM3_PER_UMOL) * V_drop)
    return {"C": C, "V_drop": V_drop}
