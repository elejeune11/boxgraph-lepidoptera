import math

# 1 uM = 1e-21 mol/um^3; N_A = 6.02214076e23 mol^-1 => 602.214076 monomers / (uM * um^3)
MONOMERS_PER_UM3_PER_UM = 602.214076


def run(inputs, settings):
    V_d = (4.0 / 3.0) * math.pi * (inputs["R_d"] ** 3)
    L_tot = inputs["c_a"] * MONOMERS_PER_UM3_PER_UM * V_d * (inputs["a_mon"] * 1e-3)
    L_B = 2.0 * math.pi * inputs["R_r"]
    n = L_tot / L_B
    return {"V_d": V_d, "L_tot": L_tot, "L_B": L_B, "n": n}


def symbolic():
    import sympy as sp

    c_a, R_d, R_r, a_mon = sp.symbols("c_a R_d R_r a_mon", positive=True)
    V_d = sp.Rational(4, 3) * sp.pi * R_d**3
    L_tot = c_a * MONOMERS_PER_UM3_PER_UM * V_d * (a_mon * 1e-3)
    L_B = 2 * sp.pi * R_r
    n = L_tot / L_B
    return {"V_d": V_d, "L_tot": L_tot, "L_B": L_B, "n": n}
