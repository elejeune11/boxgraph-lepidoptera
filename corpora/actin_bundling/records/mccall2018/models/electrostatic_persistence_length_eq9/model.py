import math

N_A = 6.02214076e23
MOLAR_TO_PER_A3 = N_A * 1e-27


def run(inputs, settings):
    n_tot = (inputs["Z_c"] ** 2 * inputs["c_c"] + inputs["S_salt"]) * MOLAR_TO_PER_A3
    return {"l_E": (inputs["f"] ** 2) / (16.0 * math.pi * n_tot * (inputs["b"] ** 2))}


def symbolic():
    import sympy as sp

    f, Z_c, c_c, S_salt, b = sp.symbols("f Z_c c_c S_salt b", positive=True)
    return {"l_E": f ** 2 / (16 * sp.pi * sp.Float(MOLAR_TO_PER_A3) * (Z_c ** 2 * c_c + S_salt) * b ** 2)}
