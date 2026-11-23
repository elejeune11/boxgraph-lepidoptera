import math


def run(inputs, settings):
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    return {"g_c": (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)}


def symbolic():
    import sympy as sp
    nu, h, L = sp.symbols("nu h L", positive=True)
    return {"g_c": sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2}
