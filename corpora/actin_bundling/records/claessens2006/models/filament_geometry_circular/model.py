import math


def run(inputs, settings):
    k_a = inputs["k_a"]
    kappa_a = inputs["kappa_a"]
    d_m = 4.0 * math.sqrt(kappa_a / k_a)
    E_Pa = (k_a * k_a) / (4.0 * math.pi * kappa_a)
    return {"d": d_m * 1e9, "E": E_Pa * 1e-9}


def symbolic():
    import sympy as sp

    k_a, kappa_a = sp.symbols("k_a kappa_a", positive=True)
    d_nm = 4 * sp.sqrt(kappa_a / k_a) * 1e9
    E_GPa = (k_a**2 / (4 * sp.pi * kappa_a)) * 1e-9
    return {"d": d_nm, "E": E_GPa}
