import math


def run(inputs, settings):
    c0 = inputs["c0"]
    rescaled = ((1.0 + c0) ** 2) / c0
    geom = 1.0 + (inputs["beta"] * math.pi / inputs["L_over_l"]) ** 2
    return {
        "Pe_crit": geom * rescaled,
        "Pe_crit_rescaled": rescaled,
    }


def symbolic():
    import sympy as sp

    c0, beta, L_over_l = sp.symbols("c0 beta L_over_l", positive=True)
    rescaled = (1 + c0) ** 2 / c0
    geom = 1 + (beta * sp.pi / L_over_l) ** 2
    return {
        "Pe_crit": geom * rescaled,
        "Pe_crit_rescaled": rescaled,
    }
