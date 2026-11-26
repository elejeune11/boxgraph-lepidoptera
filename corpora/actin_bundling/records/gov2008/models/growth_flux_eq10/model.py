import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    rho_ac = inputs["rho_ac"]
    gamma = inputs["gamma"]
    F_R = inputs["F_R"]
    beta = math.pi * (a ** 3) * rho_ac / gamma
    R_dot = beta * F_R
    N_dot = (2.0 * R / (a * a)) * R_dot
    return {"beta": beta, "R_dot": R_dot, "N_dot": N_dot}


def symbolic():
    import sympy as sp

    R, a, rho_ac, gamma, F_R = sp.symbols("R a rho_ac gamma F_R", positive=True)
    beta = sp.pi * a**3 * rho_ac / gamma
    R_dot = beta * F_R
    N_dot = (2 * R / a**2) * R_dot
    return {"beta": beta, "R_dot": R_dot, "N_dot": N_dot}
