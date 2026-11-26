"""Effective interlayer shear modulus G_shear = k_x / delta."""


def run(inputs, settings):
    return {"G_shear": inputs["k_x"] / inputs["delta"]}


def symbolic():
    import sympy as sp

    k_x, delta = sp.symbols("k_x delta", positive=True)
    return {"G_shear": k_x / delta}
