# 1 uM * 1 um^3 = 1e-6 mol/L * 1e-15 L * 6.02214076e23 mol^-1 = 602.214076 molecules
MOL_PER_UM_UM3 = 602.214076

def run(inputs, settings):
    n_actin = inputs["A_tot"] * inputs["V"] * MOL_PER_UM_UM3
    return {"N_minifil": inputs["R_ma"] * n_actin / inputs["N_mol"]}

def symbolic():
    import sympy as sp
    R_ma, A_tot, V, N_mol = sp.symbols("R_ma A_tot V N_mol", positive=True)
    return {"N_minifil": R_ma * A_tot * V * MOL_PER_UM_UM3 / N_mol}
