def run(inputs, settings):
    return {"ratio": inputs["L_c"] / inputs["l"]}


def symbolic():
    import sympy as sp
    L_c, l = sp.symbols("L_c l", positive=True)
    return {"ratio": L_c / l}
