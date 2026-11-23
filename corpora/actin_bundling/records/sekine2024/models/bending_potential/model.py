"""Angular bending potential U^bend = -k cos(theta) and restoring torque k sin(theta) (Eqs. S7, S9)."""

import math


def run(inputs, settings):
    theta = inputs["theta"]
    k_bend = inputs["k_bend"]
    return {
        "U_bend": -k_bend * math.cos(theta),
        "tau_bend": k_bend * math.sin(theta),
    }


def symbolic():
    import sympy as sp

    theta, k_bend = sp.symbols("theta k_bend", positive=True)
    return {
        "U_bend": -k_bend * sp.cos(theta),
        "tau_bend": k_bend * sp.sin(theta),
    }
