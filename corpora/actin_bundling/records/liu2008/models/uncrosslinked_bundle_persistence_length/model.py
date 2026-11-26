def run(inputs, settings):
    N = float(inputs["N"])
    L_p1 = float(inputs["L_p1"])
    return {"L_p": N * L_p1}


def symbolic():
    import sympy as sp

    N, L_p1 = sp.symbols("N L_p1", positive=True)
    return {"L_p": N * L_p1}
