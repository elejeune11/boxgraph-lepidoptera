import math


def run(inputs, settings):
    theta_0 = -inputs["t_0"] * 2.0 * math.pi / inputs["M_0"]
    theta_1 = -inputs["t_1"] * 2.0 * math.pi / inputs["M_1"]
    return {"Delta_theta": theta_1 - theta_0}


def symbolic():
    import sympy as sp

    M_0, t_0, M_1, t_1 = sp.symbols("M_0 t_0 M_1 t_1", positive=True)
    return {"Delta_theta": (-t_1 * 2 * sp.pi / M_1) - (-t_0 * 2 * sp.pi / M_0)}
