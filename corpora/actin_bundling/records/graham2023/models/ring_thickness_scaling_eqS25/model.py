import math


def run(inputs, settings):
    return {"u_m": (inputs["R"] ** 1.5) / math.sqrt(inputs["L_p"])}


def symbolic():
    import sympy as sp

    R, L_p = sp.symbols("R L_p", positive=True)
    return {"u_m": R ** sp.Rational(3, 2) / sp.sqrt(L_p)}
