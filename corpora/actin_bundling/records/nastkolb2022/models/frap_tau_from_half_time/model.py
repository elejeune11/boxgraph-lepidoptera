import math


def run(inputs, settings):
    return {"tau": inputs["t_half"] / math.log(2.0)}


def symbolic():
    import sympy as sp

    t_half = sp.symbols("t_half", positive=True)
    return {"tau": t_half / sp.log(2)}
