# Exact SI Avogadro constant (mol^-1)
N_A = 6.02214076e23
# 1 L = 1e15 um^3, and 1 M = 1e6 uM -> factor = 1e21 um^3 * uM / M
UM3_PER_L_TIMES_UM_PER_M = 1e21


def run(inputs, settings):
    # Filament length density = F*L / L_d^3 = 1 / ell_m^2 [1/um^2]
    # Monomer number density = rho_G / ell_m^2 [1/um^3]
    c_uM = (inputs["rho_G"] / (inputs["ell_m"] ** 2)) * (UM3_PER_L_TIMES_UM_PER_M / N_A)
    return {"c_G": c_uM}


def symbolic():
    import sympy as sp

    ell_m, rho_G = sp.symbols("ell_m rho_G", positive=True)
    return {"c_G": (rho_G / ell_m**2) * sp.Float(UM3_PER_L_TIMES_UM_PER_M / N_A)}
