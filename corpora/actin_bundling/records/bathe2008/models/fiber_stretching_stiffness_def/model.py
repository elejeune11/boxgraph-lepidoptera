"""Discrete fiber stretching stiffness between cross-links (p. 2 col 1): k_f := E_f * A_f / delta."""


def run(inputs, settings):
    return {"k_f": inputs["E_f_A_f"] / inputs["delta"]}


def symbolic():
    import sympy as sp

    E_f_A_f, delta = sp.symbols("E_f_A_f delta", positive=True)
    return {"k_f": E_f_A_f / delta}
