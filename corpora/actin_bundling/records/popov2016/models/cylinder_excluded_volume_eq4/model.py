import math

def run(inputs, settings):
    # Exact analytical evaluation of Eq. 4 for two aligned parallel cylinders of length l_0
    # separated by perpendicular center-to-center distance d:
    # int_0^1 int_0^1 ds dt / (d^2 + l_0^2 (t - s)^2)^2 = atan(l_0 / d) / (d^3 * l_0).
    d = inputs["d"]
    l0 = inputs["l_0"]
    return {"U_vol": inputs["K_vol"] * math.atan(l0 / d) / (d**3 * l0)}

def symbolic():
    import sympy as sp
    K_vol, l_0, d = sp.symbols("K_vol l_0 d", positive=True)
    return {"U_vol": K_vol * sp.atan(l_0 / d) / (d**3 * l_0)}
