def run(inputs, settings):
    return {"ratio": inputs["L_M"] / inputs["L_D"]}


def symbolic():
    import sympy as sp
    L_M, L_D = sp.symbols("L_M L_D", positive=True)
    return {"ratio": L_M / L_D}
