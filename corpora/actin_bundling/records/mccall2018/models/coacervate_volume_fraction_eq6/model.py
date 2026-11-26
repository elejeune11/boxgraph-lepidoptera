def run(inputs, settings):
    f_pp = inputs["f_pp"]
    rho_ratio = inputs["rho_pp"] / inputs["rho_solvent"]
    return {"Phi_pp": f_pp / (f_pp + (1.0 - f_pp) * rho_ratio)}


def symbolic():
    import sympy as sp

    f_pp, rho_pp, rho_solvent = sp.symbols("f_pp rho_pp rho_solvent", positive=True)
    return {"Phi_pp": f_pp / (f_pp + (1 - f_pp) * (rho_pp / rho_solvent))}
