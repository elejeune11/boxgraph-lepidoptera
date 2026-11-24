import math


def run(inputs, settings):
    return {"y": inputs["a"] * (1.0 - math.exp(-inputs["t"] * math.log(2.0) / inputs["t_half"])) + inputs["c"]}


def symbolic():
    import sympy as sp

    a, t_half, c, t = sp.symbols("a t_half c t", positive=True)
    return {"y": a * (1 - sp.exp(-t * sp.log(2) / t_half)) + c}
