import math


def run(inputs, settings):
    R_max = inputs["R_max"]
    R_f = inputs["R_f"]
    theta_rad = math.radians(inputs["theta"])
    p_star = 2.0 * R_f * (R_max - R_f * math.cos(theta_rad)) / (R_max**2 - R_f**2)
    return {"P_star": p_star}


def symbolic():
    import sympy as sp
    R_max, R_f, theta = sp.symbols("R_max R_f theta", positive=True)
    return {
        "P_star": 2 * R_f * (R_max - R_f * sp.cos(theta * sp.pi / 180)) / (R_max**2 - R_f**2)
    }
