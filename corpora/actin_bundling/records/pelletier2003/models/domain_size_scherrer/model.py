import math


def run(inputs, settings):
    return {"L": inputs["K"] * 2.0 * math.pi / inputs["Delta_q"]}


def symbolic():
    import sympy as sp

    Delta_q, K = sp.symbols("Delta_q K", positive=True)
    return {"L": K * 2 * sp.pi / Delta_q}
