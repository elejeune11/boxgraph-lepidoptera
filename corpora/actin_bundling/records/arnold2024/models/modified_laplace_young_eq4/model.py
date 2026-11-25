def run(inputs, settings):
    P_Ld_eff = inputs["P_Ld"] + inputs["sigma_el"]
    kappa = (inputs["P_Lo"] - P_Ld_eff) / inputs["lam"]
    return {"P_Ld_eff": P_Ld_eff, "kappa": kappa}


def symbolic():
    import sympy as sp

    P_Lo, P_Ld, sigma_el, lam = sp.symbols("P_Lo P_Ld sigma_el lam", real=True)
    P_Ld_eff = P_Ld + sigma_el
    kappa = (P_Lo - P_Ld - sigma_el) / lam
    return {"P_Ld_eff": P_Ld_eff, "kappa": kappa}
