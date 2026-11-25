def run(inputs, settings):
    return {"AR": inputs["L_D"] / inputs["L_y"]}


def symbolic():
    import sympy as sp
    L_D, L_y = sp.symbols("L_D L_y", positive=True)
    return {"AR": L_D / L_y}
