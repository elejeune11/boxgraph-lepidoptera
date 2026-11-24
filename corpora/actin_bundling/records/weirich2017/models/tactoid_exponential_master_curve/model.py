import math


def run(inputs, settings):
    return {"L_scaled": math.exp(-inputs["t_scaled"])}


def symbolic():
    import sympy as sp

    t_scaled = sp.symbols("t_scaled", real=True)
    return {"L_scaled": sp.exp(-t_scaled)}
