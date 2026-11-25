import math


def run(inputs, settings):
    R_f = inputs["R_f"]
    theta_rad = math.radians(inputs["theta"])
    c = inputs["c"]
    dp_2s = math.cos(theta_rad) / R_f + c / (R_f**2)
    return {
        "DeltaP_over_2sigma": dp_2s,
        "DeltaP_over_sigma": 2.0 * dp_2s,
    }


def symbolic():
    import sympy as sp
    R_f, theta = sp.symbols("R_f theta", positive=True)
    c = sp.Symbol("c", real=True)
    dp_2s = sp.cos(theta * sp.pi / 180) / R_f + c / R_f**2
    return {
        "DeltaP_over_2sigma": dp_2s,
        "DeltaP_over_sigma": 2 * dp_2s,
    }
