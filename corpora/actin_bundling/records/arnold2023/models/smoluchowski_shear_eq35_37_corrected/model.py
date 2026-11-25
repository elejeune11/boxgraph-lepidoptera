import math
from scipy.optimize import brentq

K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = float(inputs["T"])
    eta = float(inputs["eta"])
    phi_A = float(inputs["phi_A"])
    gamma_dot = float(inputs["gamma_dot"])
    C = float(inputs["C"])
    eps = float(inputs["epsilon"])
    a0 = float(inputs["a_0"])
    t = float(inputs["t"])
    tol = float(settings.get("root_tol", 1.0e-11))

    m = (K_B * T * phi_A / (16.0 * eta * eps)) * 1.0e18  # um^3 / s
    C1 = C * math.sqrt(16.0 * gamma_dot * eta / (K_B * T)) * 1.0e-9  # um^(-3/2)

    def F(a):
        u = C1 * (a ** 1.5)
        if u < 1.0e-5:
            return (a ** 3 / 3.0) * (1.0 - (2.0 / 3.0) * u + 0.5 * u * u)
        return (2.0 / (3.0 * C1 * C1)) * (u - math.log1p(u))

    C0 = F(a0)
    target = C0 + m * t
    hi = max(
        a0 * 2.0,
        (3.0 * target) ** (1.0 / 3.0) * 2.0,
        (1.5 * C1 * target) ** (2.0 / 3.0) * 2.0 if C1 > 0.0 else a0 * 2.0,
    )
    while F(hi) < target:
        hi *= 2.0
    a = brentq(lambda x: F(x) - target, a0, hi, xtol=tol, rtol=tol)
    u = C1 * (a ** 1.5)
    da_dt = (m / (a * a)) * (1.0 + u)
    alpha_eff = (t / a) * da_dt
    return {
        "m": float(m),
        "a": float(a),
        "da_dt": float(da_dt),
        "C1_a_three_halves": float(u),
        "alpha_eff": float(alpha_eff),
    }
