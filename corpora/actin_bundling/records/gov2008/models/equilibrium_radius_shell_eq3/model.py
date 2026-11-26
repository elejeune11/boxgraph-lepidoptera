def run(inputs, settings):
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    E_c_star = inputs["E_c_star"]
    if E_c <= E_c_star:
        return {"R_0": None}
    return {"R_0": lam / (1.0 - E_c_star / E_c)}


def symbolic():
    import sympy as sp

    lam, E_c, E_c_star = sp.symbols("lam E_c E_c_star", positive=True)
    return {"R_0": lam / (1 - E_c_star / E_c)}
