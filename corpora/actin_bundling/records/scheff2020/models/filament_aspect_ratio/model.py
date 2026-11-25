def run(inputs, settings):
    return {"a_f": inputs["L"] / inputs["D"]}


def symbolic():
    import sympy as sp

    L, D = sp.symbols("L D", positive=True)
    return {"a_f": L / D}
