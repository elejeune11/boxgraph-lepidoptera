def run(inputs, settings):
    return {"asphericity": ((inputs["c"] ** 2) / (inputs["a"] * inputs["b"])) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    a, b, c = sp.symbols("a b c", positive=True)
    return {"asphericity": (c**2 / (a * b)) ** sp.Rational(1, 3)}
