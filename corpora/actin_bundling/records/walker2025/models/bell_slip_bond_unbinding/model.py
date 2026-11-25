import math


def run(inputs, settings):
    return {
        "k_unbind": float(inputs["k_unbind_0"]) * math.exp(float(inputs["F"]) / float(inputs["F_unbind"]))
    }


def symbolic():
    import sympy as sp

    k_unbind_0, F, F_unbind = sp.symbols("k_unbind_0 F F_unbind", positive=True)
    return {"k_unbind": k_unbind_0 * sp.exp(F / F_unbind)}
