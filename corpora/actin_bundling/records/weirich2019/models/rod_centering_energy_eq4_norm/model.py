import math


def _eq4_shifted_dimless(u, r):
    a = abs(u)
    base0 = -2.0 * math.log(1.0 - 0.25 * r * r)
    if a <= 0.5 * r:
        val = -math.log(1.0 - (0.5 * r - a) ** 2) - math.log(1.0 - (0.5 * r + a) ** 2)
    else:
        val = math.log(1.0 - (a - 0.5 * r) ** 2) - math.log(1.0 - (a + 0.5 * r) ** 2)
    return val - base0


def run(inputs, settings):
    u = float(inputs["u"])
    L_D = float(inputs["L_D"])
    L_M = float(inputs["L_M"])
    r = L_M / L_D
    return {"Delta_F_shifted_over_F0": _eq4_shifted_dimless(u, r)}
