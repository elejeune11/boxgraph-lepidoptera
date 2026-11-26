"""Cross-link shear stiffness k_x = (k_x / delta) * delta."""


def run(inputs, settings):
    return {"k_x": inputs["G_shear"] * inputs["delta"]}


def symbolic():
    import sympy as sp

    G_shear, delta = sp.symbols("G_shear delta", positive=True)
    return {"k_x": G_shear * delta}
