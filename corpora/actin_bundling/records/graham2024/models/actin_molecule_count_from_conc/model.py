import math

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    return {"N": inputs["C"] * V_um3 * UM_UM3_TO_MOLECULES}


def symbolic():
    import sympy as sp

    C, R_drop = sp.symbols("C R_drop", positive=True)
    return {"N": C * sp.Rational(4, 3) * sp.pi * R_drop**3 * UM_UM3_TO_MOLECULES}
