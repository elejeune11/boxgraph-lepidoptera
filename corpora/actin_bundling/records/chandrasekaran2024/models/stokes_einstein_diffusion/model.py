import math


def run(inputs, settings):
    # kBT in pN*nm -> 1e-3 pN*um; eta in pN*s/um^2; r in nm -> 1e-3 um
    # kBT / (6*pi*eta*r) = (kBT * 1e-3) / (6*pi*eta*(r * 1e-3)) = kBT / (6*pi*eta*r) in um^2/s
    return {"D_SE": inputs["kBT"] / (6.0 * math.pi * inputs["eta"] * inputs["r"])}


def symbolic():
    import sympy as sp

    kBT, eta, r = sp.symbols("kBT eta r", positive=True)
    return {"D_SE": kBT / (6 * sp.pi * eta * r)}
