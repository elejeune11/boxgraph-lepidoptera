import math

def run(inputs, settings):
    er = inputs["eta_bar_rot"]
    bb = inputs["beta_bar"]
    return {
        "margin": 8.0 * er - bb * bb,
        "beta_bar_min": -2.0 * math.sqrt(2.0 * er),
    }

def symbolic():
    import sympy as sp
    eta_bar_rot, beta_bar = sp.symbols("eta_bar_rot beta_bar", real=True)
    return {
        "margin": 8 * eta_bar_rot - beta_bar**2,
        "beta_bar_min": -2 * sp.sqrt(2 * eta_bar_rot),
    }
