import math


def run(inputs, settings):
    return {"l_t": math.sqrt(inputs["D"] * inputs["tau"])}


def symbolic():
    import sympy as sp

    D, tau = sp.symbols("D tau", positive=True)
    return {"l_t": sp.sqrt(D * tau)}
