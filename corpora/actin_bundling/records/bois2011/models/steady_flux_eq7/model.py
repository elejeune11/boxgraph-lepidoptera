import math


def run(inputs, settings):
    rhs = inputs["delta_sigma"] / inputs["gamma"] - inputs["D"] * math.log(inputs["c_ratio"])
    return {"j0": rhs / inputs["inv_c_int"]}


def symbolic():
    import sympy as sp

    delta_sigma = sp.Symbol("delta_sigma", real=True)
    gamma, D, c_ratio, inv_c_int = sp.symbols("gamma D c_ratio inv_c_int", positive=True)
    return {"j0": (delta_sigma / gamma - D * sp.log(c_ratio)) / inv_c_int}
