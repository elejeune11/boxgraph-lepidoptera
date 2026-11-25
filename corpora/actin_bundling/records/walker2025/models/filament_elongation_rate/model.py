import math


def run(inputs, settings):
    return {"k_grow": (2.0 * math.pi * float(inputs["R_drop"]) - float(inputs["L0"])) / float(inputs["t_sim"])}


def symbolic():
    import sympy as sp

    R_drop, L0, t_sim = sp.symbols("R_drop L0 t_sim", positive=True)
    return {"k_grow": (2 * sp.pi * R_drop - L0) / t_sim}
