def run(inputs, settings):
    return {"L_B": float(inputs["L_A"])}


def symbolic():
    import sympy as sp
    L_A = sp.symbols("L_A", positive=True)
    return {"L_B": L_A}
