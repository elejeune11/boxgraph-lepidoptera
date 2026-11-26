import math


def run(inputs, settings):
    return {"C_s": math.exp(-inputs["s"] / inputs["L_p_eff"])}


def symbolic():
    import sympy as sp

    s, L_p_eff = sp.symbols("s L_p_eff", positive=True)
    return {"C_s": sp.exp(-s / L_p_eff)}
