def run(inputs, settings):
    return {"F_0": inputs["k"] * (inputs["L_D"] ** 2) * 1000.0}


def symbolic():
    import sympy as sp
    k, L_D = sp.symbols("k L_D", positive=True)
    return {"F_0": 1000 * k * L_D ** 2}
