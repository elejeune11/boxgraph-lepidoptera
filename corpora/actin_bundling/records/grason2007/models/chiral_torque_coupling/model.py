def run(inputs, settings):
    return {"gamma": 2.0 * inputs["K2"] / inputs["P"]}


def symbolic():
    import sympy as sp
    K2, P = sp.symbols("K2 P", positive=True)
    return {"gamma": 2 * K2 / P}
