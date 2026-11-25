import math


def run(inputs, settings):
    n = inputs["n"]
    theta_rad = math.radians(inputs["theta"])
    return {"P_star": 2.0 * (n - math.cos(theta_rad)) / (n**2 - 1.0)}


def symbolic():
    import sympy as sp
    n, theta = sp.symbols("n theta", positive=True)
    return {"P_star": 2 * (n - sp.cos(theta * sp.pi / 180)) / (n**2 - 1)}
