def run(inputs, settings):
    x = inputs["X_A"]
    return {"gamma": x * (x - 1.0) / 2.0}


def symbolic():
    import sympy as sp
    X_A = sp.symbols("X_A", positive=True)
    return {"gamma": X_A * (X_A - 1) / 2}
