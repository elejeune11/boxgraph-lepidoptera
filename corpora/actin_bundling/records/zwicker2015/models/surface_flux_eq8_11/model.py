import math


def run(inputs, settings):
    R = inputs["R"]
    D_B = inputs["D_B"]
    l = inputs["l"]
    gb = inputs["gamma_beta"]
    phi_0_B = inputs["phi_0_B"]
    return {"J": 4.0 * math.pi * R * D_B * (1.0 + R / l) * (gb / R - phi_0_B)}


def symbolic():
    import sympy as sp

    R, D_B, l, gamma_beta, phi_0_B = sp.symbols("R D_B l gamma_beta phi_0_B", positive=True)
    return {"J": 4 * sp.pi * R * D_B * (1 + R / l) * (gamma_beta / R - phi_0_B)}
