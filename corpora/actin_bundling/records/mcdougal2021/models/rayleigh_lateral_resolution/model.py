"""Rayleigh lateral diffraction resolution stand-in: r_lat = 0.61 * lambda_0 / NA."""


def run(inputs, settings):
    return {"r_lat": 0.61 * float(inputs["lambda_0"]) / float(inputs["NA"])}


def symbolic():
    import sympy as sp

    lambda_0, NA = sp.symbols("lambda_0 NA", positive=True)
    return {"r_lat": sp.Rational(61, 100) * lambda_0 / NA}
