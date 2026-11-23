import math


def run(inputs, settings):
    a = inputs["alpha"]
    phi = inputs["phi"]
    i = inputs["i"]
    k = i * math.pi
    num = (a * a - 2.0 * k * k) * ((k * a * a * phi + k ** 3) ** 2)
    den = 3.0 * a * a * (a * a * phi + 2.0 * k * k) * (a ** 4 * phi - 4.0 * k ** 4)
    return {"chi_i2": num / den}


def symbolic():
    import sympy as sp

    alpha, phi, i = sp.symbols("alpha phi i", positive=True)
    k = i * sp.pi
    num = (alpha**2 - 2 * k**2) * (k * alpha**2 * phi + k**3) ** 2
    den = 3 * alpha**2 * (alpha**2 * phi + 2 * k**2) * (alpha**4 * phi - 4 * k**4)
    return {"chi_i2": num / den}
