import math


def run(inputs, settings):
    return {"u_m": (inputs["R"] ** 2) / math.sqrt(inputs["L"] * inputs["L_p"])}


def symbolic():
    import sympy as sp

    R, L, L_p = sp.symbols("R L L_p", positive=True)
    return {"u_m": R ** 2 / sp.sqrt(L * L_p)}
