"""Pseudo-image cell side length L_box / N_grid (Methods p. 11)."""


def run(inputs, settings):
    return {"dx_pixel": inputs["L_box"] / inputs["N_grid"]}


def symbolic():
    import sympy as sp

    L_box, N_grid = sp.symbols("L_box N_grid", positive=True)
    return {"dx_pixel": L_box / N_grid}
