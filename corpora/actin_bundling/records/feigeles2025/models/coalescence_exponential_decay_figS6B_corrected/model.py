import math


def run(inputs, settings):
    L_i = inputs["L_i"]
    L_f = inputs["L_f"]
    tau = inputs["tau"]
    t = inputs["t"]
    return {"L": L_f + (L_i - L_f) * math.exp(-t / tau)}


def symbolic():
    import sympy as sp
    L_i, L_f, tau, t = sp.symbols("L_i L_f tau t", positive=True)
    return {"L": L_f + (L_i - L_f) * sp.exp(-t / tau)}
