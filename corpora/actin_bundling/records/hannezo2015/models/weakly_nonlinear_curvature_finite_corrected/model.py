import math


def run(inputs, settings):
    a = inputs["alpha"]
    phi = inputs["phi"]
    k = inputs["i"] * math.pi
    num = (2.0 * k * k - a * a) * (a * a + k * k) * ((a * a * phi + k * k) ** 2)
    den = 3.0 * a * a * k * k * (4.0 * k ** 4 - a ** 4 * phi)
    return {"chi_i2": num / den}


def symbolic():
    import sympy as sp

    alpha, phi, i = sp.symbols("alpha phi i", positive=True)
    k = i * sp.pi
    num = (2 * k**2 - alpha**2) * (alpha**2 + k**2) * (alpha**2 * phi + k**2) ** 2
    den = 3 * alpha**2 * k**2 * (4 * k**4 - alpha**4 * phi)
    return {"chi_i2": num / den}
