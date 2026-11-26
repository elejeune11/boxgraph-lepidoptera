def run(inputs, settings):
    return {"n_component": inputs["fraction"] * inputs["n_total"]}


def symbolic():
    import sympy as sp

    fraction, n_total = sp.symbols("fraction n_total", positive=True)
    return {"n_component": fraction * n_total}
