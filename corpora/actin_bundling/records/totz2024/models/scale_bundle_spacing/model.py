def run(inputs, settings):
    return {"L": float(inputs["m"] * inputs["L_0"])}


def symbolic():
    import sympy as sp
    m, L_0 = sp.symbols("m L_0", positive=True)
    return {"L": m * L_0}
