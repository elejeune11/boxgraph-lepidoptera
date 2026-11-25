import math


def run(inputs, settings):
    L_max = 2.0 * math.pi * inputs["R_drop"]
    return {"L_max": L_max, "N_seg": L_max / inputs["L_seg"]}


def symbolic():
    import sympy as sp

    R_drop, L_seg = sp.symbols("R_drop L_seg", positive=True)
    L_max = 2 * sp.pi * R_drop
    return {"L_max": L_max, "N_seg": L_max / L_seg}
