import math


def run(inputs, settings):
    a = inputs["a"]
    p = inputs["p"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    # Setting g = 0 in printed Eq. 3 (evaluated in SI base units: m, Pa):
    # 0 = 1 + 3*h^2*a^2/4 - p*L^2 / (a*pi^2*Y_min*g_c*h^3)
    Y_min = p * (L ** 2) / (a * (math.pi ** 2) * g_c * (h ** 3) * (1.0 + 3.0 * (h ** 2) * (a ** 2) / 4.0))
    return {"Y_min": Y_min}


def symbolic():
    import sympy as sp
    a, p, nu, h, L = sp.symbols("a p nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    return {"Y_min": p * L**2 / (a * sp.pi**2 * g_c * h**3 * (1 + 3 * h**2 * a**2 / 4))}
