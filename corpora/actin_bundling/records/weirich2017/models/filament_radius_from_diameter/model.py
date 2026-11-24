def run(inputs, settings):
    return {"r": 0.5 * inputs["d"]}


def symbolic():
    import sympy as sp

    d = sp.Symbol("d", positive=True)
    return {"r": d / 2}
