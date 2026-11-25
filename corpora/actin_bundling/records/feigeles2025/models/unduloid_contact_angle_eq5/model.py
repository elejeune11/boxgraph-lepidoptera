import math


def run(inputs, settings):
    R_f = inputs["R_f"]
    dp_2s = inputs["DeltaP_over_2sigma"]
    c = inputs["c"]
    cos_theta = max(-1.0, min(1.0, R_f * dp_2s - c / R_f))
    return {"theta": math.degrees(math.acos(cos_theta))}


def symbolic():
    import sympy as sp
    R_f, DeltaP_over_2sigma = sp.symbols("R_f DeltaP_over_2sigma", positive=True)
    c = sp.Symbol("c", real=True)
    return {"theta": (180 / sp.pi) * sp.acos(R_f * DeltaP_over_2sigma - c / R_f)}
