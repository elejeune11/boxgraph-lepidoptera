import math


def run(inputs, settings):
    K_vol, L, d = inputs["K_vol"], inputs["L"], inputs["d"]
    return {"U_vol": K_vol * (L ** 2) * math.atan(L / d) / (L * (d ** 3))}


def symbolic():
    import sympy as sp
    K_vol, L, d = sp.symbols("K_vol L d", positive=True)
    return {"U_vol": K_vol * L * sp.atan(L / d) / (d ** 3)}
