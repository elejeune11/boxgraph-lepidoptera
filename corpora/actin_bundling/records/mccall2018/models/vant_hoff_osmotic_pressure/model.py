K_B = 1.380649e-23      # J/K (exact SI)
N_A = 6.02214076e23     # 1/mol (exact SI)


def run(inputs, settings):
    # c in uM -> 1 uM = 1e-6 mol/L = 1e-3 mol/m^3
    c_mol_m3 = inputs["c"] * 1.0e-3
    return {"Pi": c_mol_m3 * N_A * K_B * inputs["T"]}


def symbolic():
    import sympy as sp

    c, T = sp.symbols("c T", positive=True)
    return {"Pi": sp.Float(1.0e-3 * N_A * K_B) * c * T}
