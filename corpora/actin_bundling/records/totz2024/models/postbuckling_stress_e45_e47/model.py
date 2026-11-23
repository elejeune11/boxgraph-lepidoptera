import math


def run(inputs, settings):
    Y = inputs["Y"]
    g = inputs["g"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    return {
        "sigma_xx": -Y * g_c,
        "sigma_yy": -Y * (g - nu * (g - g_c)),
        "sigma_xy": 0.0,
    }


def symbolic():
    import sympy as sp
    Y, g, nu, h, L = sp.symbols("Y g nu h L", positive=True)
    g_c = sp.pi**2 / (3 * (1 - nu**2)) * (h / L)**2
    return {
        "sigma_xx": -Y * g_c,
        "sigma_yy": -Y * (g - nu * (g - g_c)),
        "sigma_xy": sp.Integer(0),
    }
