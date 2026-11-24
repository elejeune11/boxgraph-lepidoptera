import math


def run(inputs, settings):
    # eta [Pa*s] * (pi * R_0^2 * lam) [um^3 = 1e-18 m^3] * (v / R_0)^2 [(1/s)^2] = W
    R_0 = inputs["R_0"]
    vol_m3 = math.pi * (R_0**2) * inputs["lam"] * 1e-18
    shear_rate = inputs["v"] / R_0
    return {"D": inputs["eta"] * vol_m3 * (shear_rate**2)}


def symbolic():
    import sympy as sp

    eta, R_0, lam, v = sp.symbols("eta R_0 lam v", positive=True)
    return {
        "D": eta
        * (sp.pi * R_0**2 * lam * sp.Rational(1, 10**18))
        * (v / R_0) ** 2
    }
