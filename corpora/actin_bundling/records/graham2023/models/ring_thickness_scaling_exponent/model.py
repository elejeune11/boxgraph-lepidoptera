import math


def _u_m(D, L_p):
    R = 0.5 * D
    L = 2.0 * math.pi * R
    return (R ** 2) / math.sqrt(L * L_p)


def run(inputs, settings):
    D = inputs["D"]
    L_p = inputs["L_p"]
    eps = 1e-6
    u1 = _u_m(D, L_p)
    u2 = _u_m(D * (1.0 + eps), L_p)
    b = math.log(u2 / u1) / math.log(1.0 + eps)
    return {"b": b}


def symbolic():
    import sympy as sp

    D, L_p = sp.symbols("D L_p", positive=True)
    R = D / 2
    L = 2 * sp.pi * R
    u = R ** 2 / sp.sqrt(L * L_p)
    return {"b": sp.simplify(D * sp.diff(sp.log(u), D))}
