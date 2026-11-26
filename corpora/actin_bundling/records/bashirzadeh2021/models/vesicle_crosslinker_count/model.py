import math


def run(inputs, settings):
    return {"N_xl": inputs["rho_xl"] * math.pi * (inputs["r_c"] ** 2)}


def symbolic():
    import sympy as sp

    rho_xl, r_c = sp.symbols("rho_xl r_c", positive=True)
    return {"N_xl": rho_xl * sp.pi * r_c ** 2}
