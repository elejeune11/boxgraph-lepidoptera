import math


def run(inputs, settings):
    return {"J": 4.0 * math.pi * inputs["D_B"] * (inputs["gamma_beta"] - inputs["R"] * inputs["phi_0_B"])}


def symbolic():
    import sympy as sp

    R, D_B, gamma_beta, phi_0_B = sp.symbols("R D_B gamma_beta phi_0_B", positive=True)
    return {"J": 4 * sp.pi * D_B * (gamma_beta - R * phi_0_B)}
