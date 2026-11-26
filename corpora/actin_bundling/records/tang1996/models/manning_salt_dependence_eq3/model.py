def run(inputs, settings):
    return {"K": inputs["c"] ** (-inputs["N"])}


def symbolic():
    import sympy as sp

    c, N = sp.symbols("c N", positive=True)
    return {"K": c ** (-N)}
