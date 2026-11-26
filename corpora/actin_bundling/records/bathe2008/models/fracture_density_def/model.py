"""Fiber fracture number density per fiber: phi = L / L_f."""


def run(inputs, settings):
    return {"phi": inputs["L"] / inputs["L_f"]}


def symbolic():
    import sympy as sp

    L, L_f = sp.symbols("L L_f", positive=True)
    return {"phi": L / L_f}
