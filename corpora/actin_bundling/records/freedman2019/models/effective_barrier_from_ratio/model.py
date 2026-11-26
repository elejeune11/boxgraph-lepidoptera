import math


def run(inputs, settings):
    return {"epsilon_predicted": -math.log(inputs["switching_ratio"])}


def symbolic():
    import sympy as sp

    switching_ratio = sp.Symbol("switching_ratio", positive=True)
    return {"epsilon_predicted": -sp.log(switching_ratio)}
