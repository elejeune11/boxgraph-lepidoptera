import math


def run(inputs, settings):
    return {
        "k_alpha_unbind": inputs["k0_alpha_unbind"]
        * math.exp(inputs["F_ext"] / inputs["F_alpha_unbind"])
    }


def symbolic():
    import sympy as sp
    k0_alpha_unbind, F_ext, F_alpha_unbind = sp.symbols(
        "k0_alpha_unbind F_ext F_alpha_unbind", positive=True
    )
    return {"k_alpha_unbind": k0_alpha_unbind * sp.exp(F_ext / F_alpha_unbind)}
