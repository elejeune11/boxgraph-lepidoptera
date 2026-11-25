import math


def run(inputs, settings):
    n = inputs["n"]
    theta_rad = math.radians(inputs["theta"])
    return {"R_FUS_over_Rf": (n**2 - 1.0) / (2.0 * (n - math.cos(theta_rad)))}


def symbolic():
    import sympy as sp
    n, theta = sp.symbols("n theta", positive=True)
    return {"R_FUS_over_Rf": (n**2 - 1) / (2 * (n - sp.cos(theta * sp.pi / 180)))}
