import math


def run(inputs, settings):
    return {"DeltaG_kBT": math.log(inputs["K_D"] / inputs["c_0"])}


def symbolic():
    import sympy as sp

    K_D, c_0 = sp.symbols("K_D c_0", positive=True)
    return {"DeltaG_kBT": sp.log(K_D / c_0)}
