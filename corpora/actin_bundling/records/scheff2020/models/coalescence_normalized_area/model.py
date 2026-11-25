import math


def run(inputs, settings):
    return {"A_norm": math.exp(-inputs["t_rescaled"])}


def symbolic():
    import sympy as sp

    t_rescaled = sp.Symbol("t_rescaled", positive=True)
    return {"A_norm": sp.exp(-t_rescaled)}
