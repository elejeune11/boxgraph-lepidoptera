def run(inputs, settings):
    return {"Pe": inputs["U"] * inputs["l"] / inputs["D"]}


def symbolic():
    import sympy as sp

    U, l, D = sp.symbols("U l D", positive=True)
    return {"Pe": U * l / D}
