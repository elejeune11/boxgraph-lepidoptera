import math


def run(inputs, settings):
    n = inputs["n"]
    ell = inputs["ell"]
    k_n = 2.0 * math.pi * n
    return {
        "k_n": k_n,
        "q_n": k_n / ell,
    }


def symbolic():
    import sympy as sp

    n, ell = sp.symbols("n ell", positive=True)
    return {
        "k_n": 2 * sp.pi * n,
        "q_n": 2 * sp.pi * n / ell,
    }
