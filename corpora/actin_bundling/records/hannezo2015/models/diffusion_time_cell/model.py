def run(inputs, settings):
    return {"T_cell": (inputs["l"] ** 2) / inputs["D"]}


def symbolic():
    import sympy as sp

    l, D = sp.symbols("l D", positive=True)
    return {"T_cell": l**2 / D}
