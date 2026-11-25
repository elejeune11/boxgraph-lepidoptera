import math


def run(inputs, settings):
    R_f = inputs["R_f"]
    R_b = inputs["R_b"]
    theta_b_rad = math.radians(inputs["theta_b"])
    theta_rad = math.radians(inputs["theta"])
    num = R_b * math.sin(theta_b_rad) - (R_b**2 * math.cos(theta_rad)) / R_f
    c = (R_f**2) * num / (R_b**2 - R_f**2)
    return {"c": c}


def symbolic():
    import sympy as sp
    R_f, R_b, theta_b, theta = sp.symbols("R_f R_b theta_b theta", positive=True)
    num = R_b * sp.sin(theta_b * sp.pi / 180) - (R_b**2 * sp.cos(theta * sp.pi / 180)) / R_f
    return {"c": R_f**2 * num / (R_b**2 - R_f**2)}
