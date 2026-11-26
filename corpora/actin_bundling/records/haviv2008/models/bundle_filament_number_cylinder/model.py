import math


def run(inputs, settings):
    D_b = float(inputs["D_b"])
    D_f = float(inputs["D_f"])
    A_b = (math.pi / 4.0) * (D_b ** 2)
    A_f = (math.pi / 4.0) * (D_f ** 2)
    N_b = A_b / A_f
    return {"A_b": A_b, "A_f": A_f, "N_b": N_b}


def symbolic():
    import sympy as sp

    D_b, D_f = sp.symbols("D_b D_f", positive=True)
    A_b = (sp.pi / 4) * D_b**2
    A_f = (sp.pi / 4) * D_f**2
    return {"A_b": A_b, "A_f": A_f, "N_b": A_b / A_f}
