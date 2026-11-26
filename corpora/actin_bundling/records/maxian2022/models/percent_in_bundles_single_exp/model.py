import math


def run(inputs, settings):
    tau1 = inputs["tau1"]
    t = inputs["t"]
    return {
        "percent": 100.0 * (1.0 - math.exp(-t / tau1)),
        "t_50": tau1 * math.log(2.0),
    }


def symbolic():
    import sympy as sp

    t, tau1 = sp.symbols("t tau1", positive=True)
    return {
        "percent": 100 * (1 - sp.exp(-t / tau1)),
        "t_50": tau1 * sp.log(2),
    }
