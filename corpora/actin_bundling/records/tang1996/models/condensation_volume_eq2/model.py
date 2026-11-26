def run(inputs, settings):
    return {"V_P": 41.1 * (inputs["zeta"] - 1.0) * (inputs["b"] ** 3)}


def symbolic():
    import sympy as sp

    zeta, b = sp.symbols("zeta b", positive=True)
    return {"V_P": sp.Rational(411, 10) * (zeta - 1) * b ** 3}
