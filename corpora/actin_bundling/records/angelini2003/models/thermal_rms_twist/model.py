import math

K_B = 1.380649e-23  # J/K (exact SI definition)


def run(inputs, settings):
    d_m = inputs["d"] * 1e-10  # angstrom -> m
    return {"Delta_theta_rms": math.sqrt(K_B * inputs["T"] * d_m / inputs["kappa"])}


def symbolic():
    import sympy as sp

    kappa, d, T = sp.symbols("kappa d T", positive=True)
    return {"Delta_theta_rms": sp.sqrt(sp.Float(K_B, 16) * T * (d * sp.Rational(1, 10**10)) / kappa)}
