import math


def run(inputs, settings):
    q_10 = inputs["q_10"]
    return {
        "q_11": math.sqrt(3.0) * q_10,
        "q_20": 2.0 * q_10,
    }


def symbolic():
    import sympy as sp

    q_10 = sp.Symbol("q_10", positive=True)
    return {
        "q_11": sp.sqrt(3) * q_10,
        "q_20": 2 * q_10,
    }
