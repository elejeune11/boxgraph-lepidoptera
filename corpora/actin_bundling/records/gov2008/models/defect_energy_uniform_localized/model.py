import math


def run(inputs, settings):
    return {"E_defect": inputs["rho"] * math.pi * (inputs["a"] ** 2) * inputs["E_c"]}


def symbolic():
    import sympy as sp

    rho, a, E_c = sp.symbols("rho a E_c", positive=True)
    return {"E_defect": rho * sp.pi * a**2 * E_c}
