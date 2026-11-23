import math


def run(inputs, settings):
    return {"n_rings": inputs["alpha"] * (inputs["phi"] ** 0.25) / (2.0 * math.pi)}


def symbolic():
    import sympy as sp

    alpha, phi = sp.symbols("alpha phi", positive=True)
    return {"n_rings": alpha * phi ** sp.Rational(1, 4) / (2 * sp.pi)}
