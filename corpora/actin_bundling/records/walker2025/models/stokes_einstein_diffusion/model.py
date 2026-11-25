import math


def run(inputs, settings):
    # eta in pN*s/um^2, r in nm (= 1e-3 um), kBT in pN*nm (= 1e-3 pN*um)
    eta = float(inputs["eta"])
    r_um = float(inputs["r"]) * 1e-3
    kBT_pN_um = float(inputs["kBT"]) * 1e-3
    mu = 1.0 / (6.0 * math.pi * eta * r_um)
    D = mu * kBT_pN_um
    return {"mu": mu, "D": D}


def symbolic():
    import sympy as sp

    eta, r, kBT = sp.symbols("eta r kBT", positive=True)
    r_um = r * sp.Rational(1, 1000)
    kBT_pN_um = kBT * sp.Rational(1, 1000)
    mu = 1 / (6 * sp.pi * eta * r_um)
    D = mu * kBT_pN_um
    return {"mu": mu, "D": D}
