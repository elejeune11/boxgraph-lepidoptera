# CODATA Avogadro constant N_A = 6.02214076e23 mol^-1
# Volume in liters: V_L = (l_comp_nm * 1e-9 m)^3 * 1e3 L/m^3 = l_comp_nm^3 * 1e-24 L
# Molecules per uM in compartment: 1e-6 mol/L * N_A * V_L = 6.02214076e-7 * l_comp_nm^3
N_A_FACTOR = 6.02214076e-7  # molecules / (uM * nm^3)

def run(inputs, settings):
    mol_per_uM = N_A_FACTOR * (inputs["l_comp"] ** 3)
    return {"k_sim": inputs["k_macro"] / mol_per_uM}

def symbolic():
    import sympy as sp
    k_macro, l_comp = sp.symbols("k_macro l_comp", positive=True)
    return {"k_sim": k_macro / (N_A_FACTOR * l_comp**3)}
