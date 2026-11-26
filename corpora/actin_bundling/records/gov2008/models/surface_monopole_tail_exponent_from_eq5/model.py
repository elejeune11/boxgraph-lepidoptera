import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_s_prime = inputs["k_s_prime"]
    rho = inputs["rho"]
    beta = inputs["beta"]
    D_R = inputs["D_R"]

    def lnP(r):
        x = r / a
        e_tot = -math.pi * x * eps_b * (x - 2.0) / (x * x) + 2.0 * math.pi * k_s_prime * (a * a) * rho * math.log(x)
        return -2.0 * beta * e_tot / D_R

    h = 1e-6
    dlnP_dlnR = (lnP(R * (1.0 + h)) - lnP(R * (1.0 - h))) / math.log((1.0 + h) / (1.0 - h))
    scale = math.pi * (a * a) * rho * beta * k_s_prime / D_R
    return {"tail_prefactor": dlnP_dlnR / scale}


def symbolic():
    import sympy as sp

    R, a, eps_b, k_s_prime, rho, beta, D_R = sp.symbols("R a eps_b k_s_prime rho beta D_R", positive=True)
    x = R / a
    e_tot = -sp.pi * x * eps_b * (x - 2) / x**2 + 2 * sp.pi * k_s_prime * a**2 * rho * sp.log(x)
    lnP = -2 * beta * e_tot / D_R
    scale = sp.pi * a**2 * rho * beta * k_s_prime / D_R
    return {"tail_prefactor": sp.simplify(R * sp.diff(lnP, R) / scale)}
