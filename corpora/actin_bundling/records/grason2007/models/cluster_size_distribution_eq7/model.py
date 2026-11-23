import math


def run(inputs, settings):
    R = inputs["R"]
    dF = inputs["dF_over_L"] * inputs["L"]
    mu = inputs["mu"]
    rho_0 = inputs["rho_0"]
    kBT = inputs["kBT"]
    exponent = -(dF - mu * rho_0 * math.pi * R**2) / kBT
    return {"weight": (R**2) * math.exp(exponent)}


def symbolic():
    import sympy as sp
    R, dF_over_L, L, mu, rho_0, kBT = sp.symbols(
        "R dF_over_L L mu rho_0 kBT", positive=True
    )
    return {
        "weight": R**2 * sp.exp(-(dF_over_L * L - mu * rho_0 * sp.pi * R**2) / kBT)
    }
