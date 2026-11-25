import math

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    L_total = inputs["N_fil"] * inputs["L_fil"]
    l_mono_um = inputs["l_mono"] * 1e-3
    N_monomers = L_total / l_mono_um
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    C_F = N_monomers / (V_um3 * UM_UM3_TO_MOLECULES)
    return {"L_total": L_total, "N_monomers": N_monomers, "C_F": C_F}


def symbolic():
    import sympy as sp

    N_fil, L_fil, l_mono, R_drop = sp.symbols("N_fil L_fil l_mono R_drop", positive=True)
    L_total = N_fil * L_fil
    N_monomers = L_total / (l_mono * sp.Rational(1, 1000))
    V_um3 = sp.Rational(4, 3) * sp.pi * R_drop**3
    C_F = N_monomers / (V_um3 * UM_UM3_TO_MOLECULES)
    return {"L_total": L_total, "N_monomers": N_monomers, "C_F": C_F}
