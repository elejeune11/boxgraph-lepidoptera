import math


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    p1 = math.sqrt(-a) if a <= 0.0 else None
    p2 = math.sqrt(-a - a_Delta) if (a + a_Delta) <= 0.0 else None
    return {"phi1_bin": p1, "phi2_bin": p2}


def symbolic():
    import sympy as sp

    a, a_Delta = sp.symbols("a a_Delta", real=True)
    return {
        "phi1_bin": sp.sqrt(-a),
        "phi2_bin": sp.sqrt(-a - a_Delta),
    }
