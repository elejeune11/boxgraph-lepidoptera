import math


def run(inputs, settings):
    u = float(inputs["u"])
    r = float(inputs["ratio"])
    a = abs(u)
    base0 = -2.0 * math.log(1.0 - 0.25 * r * r)
    if a <= 0.5 * r:
        val = -math.log(1.0 - (0.5 * r - a) ** 2) - math.log(1.0 - (0.5 * r + a) ** 2)
    else:
        val = math.log(1.0 - (a - 0.5 * r) ** 2) - math.log(1.0 - (a + 0.5 * r) ** 2)
    return {"Delta_F_shifted_over_F0": val - base0}
