import math


def run(inputs, settings):
    return {"v": inputs["k_grow"] * math.exp(-inputs["f"] / inputs["f_grow"])}


def symbolic():
    import sympy as sp

    k_grow, f, f_grow = sp.symbols("k_grow f f_grow", positive=True)
    return {"v": k_grow * sp.exp(-f / f_grow)}
