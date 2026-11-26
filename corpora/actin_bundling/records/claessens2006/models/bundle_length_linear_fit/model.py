def run(inputs, settings):
    return {"L_B": inputs["L_0"] + inputs["s"] * inputs["n"]}


def symbolic():
    import sympy as sp

    L_0, s, n = sp.symbols("L_0 s n", positive=True)
    return {"L_B": L_0 + s * n}
