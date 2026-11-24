def run(inputs, settings):
    return {"D": (inputs["a"] * (inputs["b"] ** 2)) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    a, b = sp.symbols("a b", positive=True)
    return {"D": (a * b ** 2) ** sp.Rational(1, 3)}
