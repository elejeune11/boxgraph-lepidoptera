import math


def run(inputs, settings):
    # Exact integral of K_vol * int_0^L int_0^L dx dy / (d^2 + (x-y)^2)^2
    # = K_vol * L * atan(L/d) / d^3.
    K_vol, L, d = inputs["K_vol"], inputs["L"], inputs["d"]
    return {"U_vol": K_vol * L * math.atan(L / d) / (d ** 3)}


def symbolic():
    import sympy as sp
    K_vol, L, d = sp.symbols("K_vol L d", positive=True)
    return {"U_vol": K_vol * L * sp.atan(L / d) / (d ** 3)}
