import math


def run(inputs, settings):
    k_off = inputs["k_unbind0"] * math.exp(inputs["f"] / inputs["f_u"])
    return {"k_unbind": k_off, "tau_res": 1.0 / k_off}


def symbolic():
    import sympy as sp

    k_unbind0, f, f_u = sp.symbols("k_unbind0 f f_u", positive=True)
    k_off = k_unbind0 * sp.exp(f / f_u)
    return {"k_unbind": k_off, "tau_res": 1 / k_off}
