import math


def run(inputs, settings):
    r = inputs["R_drop"]
    return {
        "L_diam": 2.0 * r,
        "L_semi": math.pi * r,
        "L_circ": 2.0 * math.pi * r,
    }


def symbolic():
    import sympy as sp

    R_drop = sp.Symbol("R_drop", positive=True)
    return {
        "L_diam": 2 * R_drop,
        "L_semi": sp.pi * R_drop,
        "L_circ": 2 * sp.pi * R_drop,
    }
