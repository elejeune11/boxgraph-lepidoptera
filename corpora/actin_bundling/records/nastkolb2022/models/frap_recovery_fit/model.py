import math


def run(inputs, settings):
    return {"y": inputs["a"] * (1.0 - math.exp(-inputs["t"] / inputs["tau"])) + inputs["c"]}


def symbolic():
    import sympy as sp

    a, tau, c, t = sp.symbols("a tau c t", positive=True)
    return {"y": a * (1 - sp.exp(-t / tau)) + c}
