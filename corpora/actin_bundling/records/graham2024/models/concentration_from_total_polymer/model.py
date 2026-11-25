import math

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    N_total = inputs["L_total"] / (inputs["l_mono"] * 1e-3)
    return {"C_T": N_total / (V_um3 * UM_UM3_TO_MOLECULES)}


def symbolic():
    import sympy as sp

    L_total, l_mono, R_drop = sp.symbols("L_total l_mono R_drop", positive=True)
    V_um3 = sp.Rational(4, 3) * sp.pi * R_drop**3
    return {"C_T": L_total / (l_mono * sp.Rational(1, 1000) * V_um3 * UM_UM3_TO_MOLECULES)}
