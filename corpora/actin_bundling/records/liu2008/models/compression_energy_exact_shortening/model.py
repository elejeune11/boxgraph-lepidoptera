import numpy as np


def run(inputs, settings):
    f_comp = float(inputs["f_comp"])    # pN
    a = float(inputs["a"])              # nm
    L = float(inputs["L"]) * 1e3        # um -> nm
    n_grid = int(settings.get("n_grid", 100))

    s = np.linspace(0.0, L, n_grid + 1)
    dr_ds = -a * (np.pi / (2.0 * L)) * np.sin(np.pi * s / (2.0 * L))
    # Numerically stable form of 1 - sqrt(1 - u) = u / (1 + sqrt(1 - u))
    u = dr_ds ** 2
    integrand = u / (1.0 + np.sqrt(np.maximum(1.0 - u, 0.0)))
    delta_L = float(np.trapezoid(integrand, s)) if hasattr(np, "trapezoid") else float(np.trapz(integrand, s))
    return {"Delta_E_comp": -f_comp * delta_L}
