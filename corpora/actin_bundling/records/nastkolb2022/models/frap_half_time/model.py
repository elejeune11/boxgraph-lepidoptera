import math


def run(inputs, settings):
    return {"t_half": inputs["tau"] * math.log(2.0)}


def symbolic():
    import sympy as sp

    tau = sp.symbols("tau", positive=True)
    return {"t_half": tau * sp.log(2)}
