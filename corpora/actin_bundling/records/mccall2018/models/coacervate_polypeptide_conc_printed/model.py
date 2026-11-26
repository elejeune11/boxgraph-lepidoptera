def run(inputs, settings):
    # rho_pp in mg/mL (= g/L); M_W in kDa (= 1000 g/mol)
    # c [mol/L] = (rho_pp * f_pp) / (1000 * M_W) -> c [mM] = (rho_pp * f_pp) / M_W
    return {"c": (inputs["rho_pp"] * inputs["f_pp"]) / inputs["M_W"]}


def symbolic():
    import sympy as sp

    rho_pp, f_pp, M_W = sp.symbols("rho_pp f_pp M_W", positive=True)
    return {"c": rho_pp * f_pp / M_W}
