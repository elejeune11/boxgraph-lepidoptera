def run(inputs, settings):
    return {"d": 2.0 * inputs["d_m"]}


def symbolic():
    import sympy as sp

    d_m = sp.Symbol("d_m", positive=True)
    return {"d": 2 * d_m}
