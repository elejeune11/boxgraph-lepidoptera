def run(inputs, settings):
    return {"b": 3.0 - 3.6 / inputs["R_over_a"]}


def symbolic():
    import sympy as sp

    R_over_a = sp.symbols("R_over_a", positive=True)
    return {"b": 3 - sp.Rational(18, 5) / R_over_a}
