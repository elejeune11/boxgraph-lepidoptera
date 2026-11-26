import math


def run(inputs, settings):
    k_a = inputs["k_a"]
    kappa_a = inputs["kappa_a"]
    d_m = math.sqrt(12.0 * kappa_a / k_a)
    E_Pa = (k_a * k_a) / (12.0 * kappa_a)
    return {"d": d_m * 1e9, "E": E_Pa * 1e-9}


def symbolic():
    import sympy as sp

    k_a, kappa_a = sp.symbols("k_a kappa_a", positive=True)
    d_nm = sp.sqrt(12 * kappa_a / k_a) * 1e9
    E_GPa = (k_a**2 / (12 * kappa_a)) * 1e-9
    return {"d": d_nm, "E": E_GPa}
