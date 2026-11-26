import math


def run(inputs, settings):
    r = inputs["r"]
    exp_eps = math.exp(inputs["eps_kBT"])
    a_um = inputs["a"] * 1e-3
    return {
        "l_F": a_um * (1.0 + r * exp_eps),
        "l_alpha": a_um * (1.0 + (1.0 / r) * exp_eps),
    }


def symbolic():
    import sympy as sp

    r, eps_kBT, a = sp.symbols("r eps_kBT a", positive=True)
    a_um = a * sp.Rational(1, 1000)
    return {
        "l_F": a_um * (1 + r * sp.exp(eps_kBT)),
        "l_alpha": a_um * (1 + sp.exp(eps_kBT) / r),
    }
