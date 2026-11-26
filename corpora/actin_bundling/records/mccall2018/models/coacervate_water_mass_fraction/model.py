def run(inputs, settings):
    return {"f_water": 100.0 * (1.0 - inputs["f_pp"])}


def symbolic():
    import sympy as sp

    f_pp = sp.Symbol("f_pp", positive=True)
    return {"f_water": 100 * (1 - f_pp)}
