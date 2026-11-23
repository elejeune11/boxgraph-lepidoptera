def run(inputs, settings):
    return {"tau_D": inputs["l"] ** 2 / inputs["D"]}


def symbolic():
    import sympy as sp

    l, D = sp.symbols("l D", positive=True)
    return {"tau_D": l**2 / D}
