import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    rho = inputs["rho"]
    beta = inputs["beta"]
    k_s_prime = inputs["k_s_prime"]
    D_R = inputs["D_R"]
    exponent = -2.0 * math.pi * (a * a) * rho * beta * k_s_prime / D_R
    return {"P": (R / a) ** exponent}


def symbolic():
    import sympy as sp

    R, a, rho, beta, k_s_prime, D_R = sp.symbols("R a rho beta k_s_prime D_R", positive=True)
    return {"P": (R / a) ** (-2 * sp.pi * a**2 * rho * beta * k_s_prime / D_R)}
