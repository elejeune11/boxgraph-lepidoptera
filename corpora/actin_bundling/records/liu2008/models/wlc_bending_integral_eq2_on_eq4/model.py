import numpy as np


def run(inputs, settings):
    kBT = float(inputs["kBT"])          # pN*nm
    L_p = float(inputs["L_p"]) * 1e3    # um -> nm
    a = float(inputs["a"])              # nm
    L = float(inputs["L"]) * 1e3        # um -> nm
    n_grid = int(settings.get("n_grid", 100))

    s = np.linspace(0.0, L, n_grid + 1)
    # Second derivative of r_x(s) = a*cos(pi*s/(2*L)) - a
    d2r_ds2 = -a * (np.pi / (2.0 * L)) ** 2 * np.cos(np.pi * s / (2.0 * L))
    integral = float(np.trapezoid(d2r_ds2 ** 2, s)) if hasattr(np, "trapezoid") else float(np.trapz(d2r_ds2 ** 2, s))
    E_bend = kBT * (L_p / 2.0) * integral
    return {"E_bend": E_bend}
