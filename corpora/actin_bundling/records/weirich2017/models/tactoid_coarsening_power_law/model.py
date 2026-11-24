def run(inputs, settings):
    return {"L": inputs["t"] ** inputs["alpha"]}


def symbolic():
    import sympy as sp

    t, alpha = sp.symbols("t alpha", positive=True)
    return {"L": t**alpha}
