import math


def run(inputs, settings):
    val = inputs["abs_ln_C0"]
    return {"C_0_upper": math.exp(val), "C_0_lower": math.exp(-val)}


def symbolic():
    import sympy as sp

    abs_ln_C0 = sp.symbols("abs_ln_C0", nonnegative=True)
    return {"C_0_upper": sp.exp(abs_ln_C0), "C_0_lower": sp.exp(-abs_ln_C0)}
