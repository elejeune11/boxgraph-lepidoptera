import math


def run(inputs, settings):
    g = inputs["g"]
    a = inputs["a"]
    p = inputs["p"]
    Y = inputs["Y"]
    nu = inputs["nu"]
    h = inputs["h"]
    # Dimensionally corrected Eq. 3 / E49 inverted for L^2:
    # A * L^4 + g * L^2 - C = 0
    A = p / (a * (math.pi ** 2) * Y * h)
    C = (1.0 + 0.75 * ((a / h) ** 2)) * (math.pi ** 2) * (h ** 2) / (3.0 * (1.0 - nu ** 2))
    L2 = (2.0 * C) / (g + math.sqrt(g * g + 4.0 * A * C))
    return {"L": float(math.sqrt(L2))}


def symbolic():
    import sympy as sp
    g, a, p, Y, nu, h = sp.symbols("g a p Y nu h", positive=True)
    A = p / (a * sp.pi**2 * Y * h)
    C = (1 + sp.Rational(3, 4) * (a / h)**2) * sp.pi**2 * h**2 / (3 * (1 - nu**2))
    return {"L": sp.sqrt(2 * C / (g + sp.sqrt(g**2 + 4 * A * C)))}
