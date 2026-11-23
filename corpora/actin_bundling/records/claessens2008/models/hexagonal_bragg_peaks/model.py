import math


def run(inputs, settings):
    q10 = inputs["q_10"]
    return {
        "q_11": math.sqrt(3.0) * q10,
        "q_20": 2.0 * q10,
        "q_21": math.sqrt(7.0) * q10,
    }


def symbolic():
    import sympy as sp

    q_10 = sp.Symbol("q_10", positive=True)
    return {
        "q_11": sp.sqrt(3) * q_10,
        "q_20": 2 * q_10,
        "q_21": sp.sqrt(7) * q_10,
    }
