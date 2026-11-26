import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    K = inputs["K"]
    rho = inputs["rho"]
    L_p = inputs["L_p"]
    beta = inputs["beta"]
    D_R = inputs["D_R"]
    return {"P": math.exp(-2.0 * beta * K * rho * (a * a) * (R * R) / (D_R * (L_p ** 4)))}


def symbolic():
    import sympy as sp

    R, a, K, rho, L_p, beta, D_R = sp.symbols("R a K rho L_p beta D_R", positive=True)
    return {"P": sp.exp(-2 * beta * K * rho * a**2 * R**2 / (D_R * L_p**4))}
