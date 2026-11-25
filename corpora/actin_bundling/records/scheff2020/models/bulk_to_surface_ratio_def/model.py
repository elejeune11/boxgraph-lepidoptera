def run(inputs, settings):
    return {"c": inputs["K_over_gamma"] / (inputs["V"] ** (1.0 / 3.0))}


def symbolic():
    import sympy as sp

    K_over_gamma, V = sp.symbols("K_over_gamma V", positive=True)
    return {"c": K_over_gamma / V ** sp.Rational(1, 3)}
