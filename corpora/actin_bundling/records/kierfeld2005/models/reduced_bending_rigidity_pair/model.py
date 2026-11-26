def run(inputs, settings):
    k1, k2 = inputs["kappa_1"], inputs["kappa_2"]
    return {"kappa_r_12": k1 * k2 / (k1 + k2)}


def symbolic():
    import sympy as sp
    kappa_1, kappa_2 = sp.symbols("kappa_1 kappa_2", positive=True)
    return {"kappa_r_12": kappa_1 * kappa_2 / (kappa_1 + kappa_2)}
