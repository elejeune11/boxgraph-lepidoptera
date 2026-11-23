import math


def run(inputs, settings):
    x = inputs["x"]
    g = inputs["g"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    if g < g_c:
        return {"w": 0.0}
    w = (h / math.sqrt(3.0)) * math.sqrt(g / g_c - 1.0) * (1.0 + math.cos(2.0 * math.pi * x / L))
    return {"w": w}


def symbolic():
    import sympy as sp
    x, g, nu, h, L = sp.symbols("x g nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    expr = sp.Piecewise(
        ((h / sp.sqrt(3)) * sp.sqrt(g / g_c - 1) * (1 + sp.cos(2 * sp.pi * x / L)), g >= g_c),
        (0, True),
    )
    return {"w": expr}
