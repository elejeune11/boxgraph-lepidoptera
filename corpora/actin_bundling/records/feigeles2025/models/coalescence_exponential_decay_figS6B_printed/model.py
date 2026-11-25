import math


def run(inputs, settings):
    # Verbatim as printed inside Fig. S6B: L = L_i + (L_i - L_f) * exp(-tau / t)
    L_i = inputs["L_i"]
    L_f = inputs["L_f"]
    tau = inputs["tau"]
    t = inputs["t"]
    return {"L": L_i + (L_i - L_f) * math.exp(-tau / t)}


def symbolic():
    import sympy as sp
    L_i, L_f, tau, t = sp.symbols("L_i L_f tau t", positive=True)
    return {"L": L_i + (L_i - L_f) * sp.exp(-tau / t)}
