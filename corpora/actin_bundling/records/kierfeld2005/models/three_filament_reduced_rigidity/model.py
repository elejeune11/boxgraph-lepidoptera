def run(inputs, settings):
    k1, k2, k3 = inputs["kappa_1"], inputs["kappa_2"], inputs["kappa_3"]
    kc12 = k1 + k2
    return {"kappa_123": kc12 * k3 / (kc12 + k3)}


def symbolic():
    import sympy as sp
    kappa_1, kappa_2, kappa_3 = sp.symbols("kappa_1 kappa_2 kappa_3", positive=True)
    kc12 = kappa_1 + kappa_2
    return {"kappa_123": kc12 * kappa_3 / (kc12 + kappa_3)}
