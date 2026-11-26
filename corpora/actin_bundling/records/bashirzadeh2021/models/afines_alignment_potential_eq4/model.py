import math


def run(inputs, settings):
    return {"U_xl_align": inputs["k_xl_align"] * (1.0 - math.cos(inputs["theta"]))}


def symbolic():
    import sympy as sp

    k_xl_align, theta = sp.symbols("k_xl_align theta", real=True)
    return {"U_xl_align": k_xl_align * (1 - sp.cos(theta))}
