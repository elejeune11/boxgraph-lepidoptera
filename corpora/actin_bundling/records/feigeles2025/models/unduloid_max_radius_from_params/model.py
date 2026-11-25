import math


def run(inputs, settings):
    R_f = inputs["R_f"]
    a = inputs["DeltaP_over_2sigma"]
    c = inputs["c"]
    disc = 1.0 + 4.0 * a * c
    R_max = (1.0 + math.sqrt(disc)) / (2.0 * a)
    return {
        "R_max": R_max,
        "n": R_max / R_f,
    }


def symbolic():
    import sympy as sp
    R_f, DeltaP_over_2sigma = sp.symbols("R_f DeltaP_over_2sigma", positive=True)
    c = sp.Symbol("c", real=True)
    R_max = (1 + sp.sqrt(1 + 4 * DeltaP_over_2sigma * c)) / (2 * DeltaP_over_2sigma)
    return {
        "R_max": R_max,
        "n": R_max / R_f,
    }
