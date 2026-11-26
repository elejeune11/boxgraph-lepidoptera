import math


def run(inputs, settings):
    return {"U_boundary": inputs["eps_boundary"] * math.exp(-inputs["d"] / inputs["lam"])}


def symbolic():
    import sympy as sp
    eps_boundary, d, lam = sp.symbols("eps_boundary d lam", positive=True)
    return {"U_boundary": eps_boundary * sp.exp(-d / lam)}
