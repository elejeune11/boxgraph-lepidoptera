import math


def run(inputs, settings):
    return {"C_s": inputs["C_0"] * math.exp(-inputs["s"] / inputs["L_p_eff"])}


def symbolic():
    import sympy as sp

    C_0, s, L_p_eff = sp.symbols("C_0 s L_p_eff", positive=True)
    return {"C_s": C_0 * sp.exp(-s / L_p_eff)}
