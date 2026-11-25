"""Phase wrapping operator W (Eqs. 1-2, McDougal et al. 2021, p. 7 & SI p. 11).

W phi(n) = phi(n) + 2*pi*k(n), with integer k(n) such that -pi < W phi(n) <= pi.
"""

import math


def wrap_phase_scalar(phi: float) -> tuple[float, int]:
    """Return (W(phi), k) such that -pi < phi + 2*pi*k <= pi."""
    k = math.floor((math.pi - phi) / (2.0 * math.pi))
    wrapped = phi + 2.0 * math.pi * k
    # Guard floating-point edge at -pi so half-open (-pi, pi] holds strictly.
    if wrapped <= -math.pi:
        k += 1
        wrapped = phi + 2.0 * math.pi * k
    elif wrapped > math.pi:
        k -= 1
        wrapped = phi + 2.0 * math.pi * k
    return wrapped, int(k)


def run(inputs, settings):
    phi = float(inputs["phi_unwrapped"])
    wrapped, k = wrap_phase_scalar(phi)
    return {"phi_wrapped": wrapped, "k": float(k)}


def symbolic():
    import sympy as sp

    phi_unwrapped = sp.Symbol("phi_unwrapped", real=True)
    k_expr = sp.floor((sp.pi - phi_unwrapped) / (2 * sp.pi))
    return {
        "phi_wrapped": phi_unwrapped + 2 * sp.pi * k_expr,
        "k": k_expr,
    }
