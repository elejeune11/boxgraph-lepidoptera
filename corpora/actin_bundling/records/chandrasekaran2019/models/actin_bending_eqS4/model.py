import math


def run(inputs, settings):
    return {"U_bend": inputs["eps_bend"] * (1.0 - math.cos(inputs["theta"]))}


def symbolic():
    import sympy as sp
    eps_bend, theta = sp.symbols("eps_bend theta", real=True)
    return {"U_bend": eps_bend * (1 - sp.cos(theta))}
