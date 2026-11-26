"""Eq. 1: alpha := k_x * L^2 / (k_f * delta^2)."""


def run(inputs, settings):
    k_x = inputs["k_x"]
    L = inputs["L"]
    k_f = inputs["k_f"]
    delta = inputs["delta"]
    return {"alpha": k_x * (L ** 2) / (k_f * (delta ** 2))}


def symbolic():
    import sympy as sp

    k_x, L, k_f, delta = sp.symbols("k_x L k_f delta", positive=True)
    return {"alpha": k_x * L**2 / (k_f * delta**2)}
