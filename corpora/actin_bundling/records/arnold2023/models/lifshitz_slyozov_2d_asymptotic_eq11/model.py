def run(inputs, settings):
    sigma = inputs["sigma"]
    t = inputs["t"]
    a = (4.5 * sigma * t) ** (1.0 / 3.0)
    return {"a": float(a), "alpha": 1.0 / 3.0}


def symbolic():
    import sympy as sp
    sigma, t = sp.symbols("sigma t", positive=True)
    return {
        "a": (sp.Rational(9, 2) * sigma * t) ** sp.Rational(1, 3),
        "alpha": sp.Rational(1, 3),
    }
