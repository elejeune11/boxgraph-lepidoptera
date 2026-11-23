"""Square ROI physical side length L_roi = N_px * dx_pixel (Methods p. 10 & Supp. Movies 2-3)."""


def run(inputs, settings):
    return {"L_roi": inputs["N_px"] * inputs["dx_pixel"]}


def symbolic():
    import sympy as sp

    N_px, dx_pixel = sp.symbols("N_px dx_pixel", positive=True)
    return {"L_roi": N_px * dx_pixel}
