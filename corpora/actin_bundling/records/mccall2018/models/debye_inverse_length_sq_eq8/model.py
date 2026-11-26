import math

N_A = 6.02214076e23          # 1/mol (exact SI)
MOLAR_TO_PER_A3 = N_A * 1e-27  # 1 mol/L = 10^-27 mol/A^3 = N_A * 10^-27 ions/A^3


def run(inputs, settings):
    n_tot = (inputs["Z_c"] ** 2 * inputs["c_c"] + inputs["S_salt"]) * MOLAR_TO_PER_A3
    return {"kappa_sq": 4.0 * math.pi * inputs["l_B"] * n_tot}


def symbolic():
    import sympy as sp

    l_B, Z_c, c_c, S_salt = sp.symbols("l_B Z_c c_c S_salt", positive=True)
    return {"kappa_sq": 4 * sp.pi * sp.Float(MOLAR_TO_PER_A3) * l_B * (Z_c ** 2 * c_c + S_salt)}
