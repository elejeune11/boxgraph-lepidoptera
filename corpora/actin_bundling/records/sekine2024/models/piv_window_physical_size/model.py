"""PIV interrogation window physical size (W_px / M_enlarge) * dx_pixel (Methods p. 10)."""


def run(inputs, settings):
    return {"W_phys": (inputs["W_px"] / inputs["M_enlarge"]) * inputs["dx_pixel"]}


def symbolic():
    import sympy as sp

    W_px, M_enlarge, dx_pixel = sp.symbols("W_px M_enlarge dx_pixel", positive=True)
    return {"W_phys": (W_px / M_enlarge) * dx_pixel}
