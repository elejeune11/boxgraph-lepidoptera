import math


def run(inputs, settings):
    phi_0_B = inputs["phi_0_B"]
    phi_plus_B = inputs["phi_plus_B"]
    R = inputs["R"]
    r = inputs["r"]
    l = inputs["l"]
    return {"phi_B": phi_0_B + (phi_plus_B - phi_0_B) * (R / r) * math.exp((R - r) / l)}


def symbolic():
    import sympy as sp

    phi_0_B, phi_plus_B, R, r, l = sp.symbols("phi_0_B phi_plus_B R r l", positive=True)
    return {"phi_B": phi_0_B + (phi_plus_B - phi_0_B) * (R / r) * sp.exp((R - r) / l)}
