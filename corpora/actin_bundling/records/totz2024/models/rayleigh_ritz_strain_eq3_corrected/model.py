import math


def run(inputs, settings):
    a = inputs["a"]
    p = inputs["p"]
    Y = inputs["Y"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    g = g_c * (1.0 + 3.0 * (a ** 2) / (4.0 * (h ** 2))) - p * (L ** 2) / (a * (math.pi ** 2) * Y * h)
    return {"g": g}


def symbolic():
    import sympy as sp
    a, p, Y, nu, h, L = sp.symbols("a p Y nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    return {"g": g_c * (1 + 3 * a**2 / (4 * h**2) - p * L**2 / (a * sp.pi**2 * Y * g_c * h))}
