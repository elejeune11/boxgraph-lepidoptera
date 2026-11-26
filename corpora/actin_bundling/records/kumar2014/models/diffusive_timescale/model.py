def run(inputs, settings):
    return {"tau": inputs["ell"] ** 2 / inputs["D"]}


def symbolic():
    import sympy as sp

    ell, D = sp.symbols("ell D", positive=True)
    return {"tau": ell**2 / D}
