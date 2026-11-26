import math

N_A = 6.02214076e23
MOLAR_TO_PER_A3 = N_A * 1e-27


def run(inputs, settings):
    A_E = 1.0 / (16.0 * math.pi * MOLAR_TO_PER_A3 * (inputs["b"] ** 2))
    return {"A_E": A_E, "l_E": A_E / inputs["c_m"]}


def symbolic():
    import sympy as sp

    c_m, b = sp.symbols("c_m b", positive=True)
    A_E = 1 / (16 * sp.pi * sp.Float(MOLAR_TO_PER_A3) * b ** 2)
    return {"A_E": A_E, "l_E": A_E / c_m}
