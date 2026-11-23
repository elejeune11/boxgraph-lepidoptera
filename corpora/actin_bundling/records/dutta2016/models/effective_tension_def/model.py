def run(inputs, settings):
    return {"epsilon": 3.0 * inputs["kappa_0"]}


def symbolic():
    import sympy as sp
    kappa_0 = sp.Symbol("kappa_0", positive=True)
    return {"epsilon": 3 * kappa_0}
