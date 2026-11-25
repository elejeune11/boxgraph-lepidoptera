def run(inputs, settings):
    return {"sigma_bar_norm": float(inputs["rho_rel"])}


def symbolic():
    import sympy as sp
    rho_rel = sp.Symbol("rho_rel", positive=True)
    return {"sigma_bar_norm": rho_rel}
