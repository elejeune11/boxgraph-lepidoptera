import math


def run(inputs, settings):
    g = inputs["g"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    if g < g_c:
        return {"a": 0.0}
    return {"a": (2.0 * h / math.sqrt(3.0)) * math.sqrt(g / g_c - 1.0)}


def symbolic():
    import sympy as sp
    g, nu, h, L = sp.symbols("g nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    expr = sp.Piecewise(
        ((2 * h / sp.sqrt(3)) * sp.sqrt(g / g_c - 1), g >= g_c),
        (0, True),
    )
    return {"a": expr}
