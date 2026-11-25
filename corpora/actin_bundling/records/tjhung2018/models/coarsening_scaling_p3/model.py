def run(inputs, settings):
    t = float(inputs["t"])
    return {"R": t ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    t = sp.Symbol("t", positive=True)
    return {"R": t ** sp.Rational(1, 3)}
