import math
from scipy.optimize import brentq


def run(inputs, settings):
    x = float(inputs["x"])          # nm
    L = float(inputs["L"])          # nm
    D = float(inputs["D"])          # nm
    L_0 = float(inputs["L_0"])      # nm
    L_p = float(inputs["L_p"])      # um
    kBT = float(inputs["kBT"])      # pN*nm
    rtol = float(settings.get("rtol", 1e-12))

    x_max = 2.0 * L / math.pi
    if x < 0.0 or x > x_max * (1.0 + 1e-12):
        raise ValueError(f"Lateral deflection x={x} nm outside geometric range [0, 2L/pi={x_max:.6g} nm]")

    if x == 0.0:
        theta = 0.0
        y = 0.0
    elif abs(x - x_max) <= 1e-12 * L:
        theta = math.pi / 2.0
        y = L * (1.0 - 2.0 / math.pi)
    else:
        target = x / L
        theta = brentq(lambda th: (1.0 - math.cos(th)) / th - target, 1e-12, math.pi / 2.0, rtol=rtol)
        y = L * (1.0 - math.sin(theta) / theta)

    d = D - 2.0 * x
    L_prot = L - L_0 - y
    E_fil = kBT * (L_p * 1e3 / (2.0 * L)) * (theta ** 2)
    return {
        "theta": theta,
        "y": y,
        "d": d,
        "L_prot": L_prot,
        "E_fil": E_fil,
    }
