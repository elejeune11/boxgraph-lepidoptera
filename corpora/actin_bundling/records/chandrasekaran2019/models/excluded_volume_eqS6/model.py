import math


def run(inputs, settings):
    # Exact integral of K_vol * int_0^1 int_0^1 ds dt / (d^2 + L^2*(t-s)^2)^2
    # = K_vol * atan(L/d) / (L * d^3).
    K_vol, L, d = inputs["K_vol"], inputs["L"], inputs["d"]
    return {"U_vol": K_vol * math.atan(L / d) / (L * (d ** 3))}


def symbolic():
    import sympy as sp
    K_vol, L, d = sp.symbols("K_vol L d", positive=True)
    return {"U_vol": K_vol * sp.atan(L / d) / (L * (d ** 3))}
