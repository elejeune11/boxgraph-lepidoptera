import math
from scipy.optimize import brentq

K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = float(inputs["T"])
    eta = float(inputs["eta"])
    phi_A = float(inputs["phi_A"])
    L = float(inputs["L"])
    a0 = float(inputs["a_0"])
    t = float(inputs["t"])
    tol = float(settings.get("root_tol", 1.0e-10))

    # Rate coefficient in um^3 / s: (m^3 -> 1e18 um^3)
    rate_um3s = (K_B * T * phi_A / (16.0 * eta)) * 1.0e18

    def G(a):
        return (a ** 3 / 9.0) * (3.0 * math.log(L / (2.0 * a)) + 1.0)

    target = G(a0) + rate_um3s * t
    a_max = 0.5 * L
    max_val = (a_max ** 3) / 9.0
    if target >= max_val:
        a = a_max * (1.0 - 1.0e-9)
    else:
        a = brentq(lambda x: G(x) - target, a0, a_max * (1.0 - 1.0e-12), xtol=tol, rtol=tol)

    da_dt = rate_um3s / (a * a * math.log(L / (2.0 * a)))
    alpha_eff = (t / a) * da_dt
    return {
        "a": float(a),
        "da_dt": float(da_dt),
        "alpha_eff": float(alpha_eff),
    }
