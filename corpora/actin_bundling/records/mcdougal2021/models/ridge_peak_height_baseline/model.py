"""Ridge peak height above the line connecting adjacent troughs (p. 7, SI p. 11, Fig. S10C).

x coordinates are in um (converted to nm via 1000 for physical Euclidean slope);
h coordinates are in nm.
"""

import math


def run(inputs, settings):
    xl = float(inputs["x_left"]) * 1000.0
    hl = float(inputs["h_left"])
    xp = float(inputs["x_peak"]) * 1000.0
    hp = float(inputs["h_peak"])
    xr = float(inputs["x_right"]) * 1000.0
    hr = float(inputs["h_right"])

    dx = xr - xl
    dh = hr - hl
    h_base = hl + ((xp - xl) / dx) * dh
    h_vert = hp - h_base
    slope = dh / dx
    h_perp = abs(h_vert) / math.sqrt(1.0 + slope * slope)
    return {"H_vert": h_vert, "H_perp": h_perp}


def symbolic():
    import sympy as sp

    x_left, h_left, x_peak, h_peak, x_right, h_right = sp.symbols(
        "x_left h_left x_peak h_peak x_right h_right", real=True
    )
    dx_nm = 1000 * (x_right - x_left)
    dh_nm = h_right - h_left
    h_vert = h_peak - (h_left + (x_peak - x_left) / (x_right - x_left) * dh_nm)
    h_perp = sp.Abs(h_vert) / sp.sqrt(1 + (dh_nm / dx_nm) ** 2)
    return {"H_vert": h_vert, "H_perp": h_perp}
