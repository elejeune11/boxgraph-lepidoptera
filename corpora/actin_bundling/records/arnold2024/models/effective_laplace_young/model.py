def run(inputs, settings):
    kappa = (inputs["P_Lo"] - inputs["P_Ld_eff"]) / inputs["lam"]
    return {"kappa": kappa}


def symbolic():
    import sympy as sp

    P_Lo, P_Ld_eff, lam = sp.symbols("P_Lo P_Ld_eff lam", real=True)
    return {"kappa": (P_Lo - P_Ld_eff) / lam}
