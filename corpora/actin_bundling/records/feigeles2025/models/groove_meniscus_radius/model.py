import math


def run(inputs, settings):
    R_max = inputs["R_max"]
    R_f = inputs["R_f"]
    theta_rad = math.radians(inputs["theta"])
    r_fus = (R_max**2 - R_f**2) / (2.0 * (R_max - R_f * math.cos(theta_rad)))
    return {
        "R_FUS": r_fus,
        "R_FUS_over_Rf": r_fus / R_f,
    }


def symbolic():
    import sympy as sp
    R_max, R_f, theta = sp.symbols("R_max R_f theta", positive=True)
    r_fus = (R_max**2 - R_f**2) / (2 * (R_max - R_f * sp.cos(theta * sp.pi / 180)))
    return {
        "R_FUS": r_fus,
        "R_FUS_over_Rf": r_fus / R_f,
    }
