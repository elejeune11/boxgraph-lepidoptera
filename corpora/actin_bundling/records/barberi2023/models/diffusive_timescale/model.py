def run(inputs, settings):
    lam = float(inputs["lam"])
    D_i = float(inputs["D_i"])
    return {"tau": (lam ** 2) / D_i}


def symbolic():
    import sympy as sp

    lam, D_i = sp.symbols("lam D_i", positive=True)
    return {"tau": (lam ** 2) / D_i}
