"""Local wrapped phase difference W Delta phi_wrapped(n) (Eq. 4, p. 7 & SI p. 11)."""

import math


def _wrap(phi: float) -> float:
    k = math.floor((math.pi - phi) / (2.0 * math.pi))
    w = phi + 2.0 * math.pi * k
    if w <= -math.pi:
        w += 2.0 * math.pi
    elif w > math.pi:
        w -= 2.0 * math.pi
    return w


def run(inputs, settings):
    diff = float(inputs["phi_n"]) - float(inputs["phi_prev"])
    return {"dphi_wrapped": _wrap(diff)}


def symbolic():
    import sympy as sp

    phi_n, phi_prev = sp.symbols("phi_n phi_prev", real=True)
    diff = phi_n - phi_prev
    k_expr = sp.floor((sp.pi - diff) / (2 * sp.pi))
    return {"dphi_wrapped": diff + 2 * sp.pi * k_expr}
