import math


def run(inputs, settings):
    return {"Delta_F": inputs["F_0"] * abs(math.sin(inputs["alpha"]))}


def symbolic():
    import sympy as sp
    F_0, alpha = sp.symbols("F_0 alpha", positive=True)
    return {"Delta_F": F_0 * sp.Abs(sp.sin(alpha))}
