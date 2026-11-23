"""Nanocluster motion anisotropy ratio <|Delta y|> / <|Delta x|> (Fig. 4d, 5a,d & Table 5)."""


def run(inputs, settings):
    return {"R_aniso": inputs["mean_abs_dy"] / inputs["mean_abs_dx"]}


def symbolic():
    import sympy as sp

    mean_abs_dx, mean_abs_dy = sp.symbols("mean_abs_dx mean_abs_dy", positive=True)
    return {"R_aniso": mean_abs_dy / mean_abs_dx}
