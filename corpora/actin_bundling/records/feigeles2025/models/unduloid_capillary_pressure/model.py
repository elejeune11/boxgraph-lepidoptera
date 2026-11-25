import math


def run(inputs, settings):
    # (uN/m) / um = (1e-6 N/m) / (1e-6 m) = 1 N/m^2
    sigma = inputs["sigma"]
    R_max = inputs["R_max"]
    R_f = inputs["R_f"]
    theta_rad = math.radians(inputs["theta"])
    dp = 2.0 * sigma * (R_max - R_f * math.cos(theta_rad)) / (R_max**2 - R_f**2)
    return {"DeltaP_condensate": dp}


def symbolic():
    import sympy as sp
    sigma, R_max, R_f, theta = sp.symbols("sigma R_max R_f theta", positive=True)
    return {
        "DeltaP_condensate": 2 * sigma * (R_max - R_f * sp.cos(theta * sp.pi / 180)) / (R_max**2 - R_f**2)
    }
