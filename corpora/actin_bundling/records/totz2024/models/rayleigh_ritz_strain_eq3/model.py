import math


def run(inputs, settings):
    a = inputs["a"]
    p = inputs["p"]
    Y = inputs["Y"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    # Evaluated literally as printed in Eq. 3 / Eq. E49 in SI base units (m, Pa):
    # g / g_c = 1 + 3*h^2*a^2/4 - p*L^2 / (a*pi^2*Y*g_c*h^3)
    g = g_c * (1.0 + 3.0 * (h ** 2) * (a ** 2) / 4.0) - p * (L ** 2) / (a * (math.pi ** 2) * Y * (h ** 3))
    return {"g": g}


def symbolic():
    import sympy as sp
    a, p, Y, nu, h, L = sp.symbols("a p Y nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    return {"g": g_c * (1 + 3 * h**2 * a**2 / 4 - p * L**2 / (a * sp.pi**2 * Y * g_c * h**3))}
