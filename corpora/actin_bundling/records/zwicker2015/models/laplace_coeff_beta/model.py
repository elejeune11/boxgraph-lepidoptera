def run(inputs, settings):
    return {"beta": 2.0 / (inputs["b"] * inputs["psi"])}


def symbolic():
    import sympy as sp

    b, psi = sp.symbols("b psi", positive=True)
    return {"beta": 2 / (b * psi)}
