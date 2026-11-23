def run(inputs, settings):
    return {"l": inputs["d"]}


def symbolic():
    import sympy as sp
    d = sp.Symbol("d", positive=True)
    return {"l": d}
