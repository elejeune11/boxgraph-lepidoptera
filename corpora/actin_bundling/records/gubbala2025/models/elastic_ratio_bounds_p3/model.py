def run(inputs, settings):
    S = float(inputs["S"])
    return {
        "E3_over_E1_min": -1.5 * S,
        "E3_over_E1_max": 3.0 / S,
    }


def symbolic():
    import sympy as sp

    S = sp.symbols("S", positive=True)
    return {
        "E3_over_E1_min": -sp.Rational(3, 2) * S,
        "E3_over_E1_max": 3 / S,
    }
