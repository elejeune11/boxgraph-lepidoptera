"""2D wrapped phase gradient magnitude and slope orientation (p. 7, SI p. 11, SI Fig. S2C-F)."""

import math


def run(inputs, settings):
    gx = float(inputs["g_x"])
    gy = float(inputs["g_y"])
    return {
        "grad_mag": math.hypot(gx, gy),
        "grad_dir": math.atan2(gy, gx),
    }


def symbolic():
    import sympy as sp

    g_x, g_y = sp.symbols("g_x g_y", real=True)
    return {
        "grad_mag": sp.sqrt(g_x**2 + g_y**2),
        "grad_dir": sp.atan2(g_y, g_x),
    }
