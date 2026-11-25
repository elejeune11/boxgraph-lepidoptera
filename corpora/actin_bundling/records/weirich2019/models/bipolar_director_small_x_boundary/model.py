import math


def run(inputs, settings):
    return {"sin_alpha_d": (inputs["x"] / inputs["L_D"]) * math.sin(inputs["alpha"])}


def symbolic():
    import sympy as sp
    x, L_D, alpha = sp.symbols("x L_D alpha", positive=True)
    return {"sin_alpha_d": (x / L_D) * sp.sin(alpha)}
