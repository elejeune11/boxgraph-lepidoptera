import math


def run(inputs, settings):
    f_m = 2.0 + inputs["m"] ** 4
    return {"u_m": (inputs["R"] ** 2) / math.sqrt(f_m * inputs["L"] * inputs["L_p"])}


def symbolic():
    import sympy as sp

    R, L, L_p, m = sp.symbols("R L L_p m", positive=True)
    return {"u_m": R ** 2 / sp.sqrt((2 + m ** 4) * L * L_p)}
