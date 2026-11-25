import math
from scipy.optimize import brentq

K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = float(inputs["T"])
    eta = float(inputs["eta"])
    gamma_dot = float(inputs["gamma_dot"])
    C = float(inputs["C"])
    eps = float(inputs["epsilon"])
    a0 = float(inputs["a_0"])
    t = float(inputs["t"])
    tol = float(settings.get("root_tol", 1.0e-11))

    # Printed Eq. (36): omits phi_A/(2*pi) and uses ln(1/2) + eps < 0 -> m < 0!
    m = (math.pi * K_B * T / (8.0 * eta * (math.log(0.5) + eps))) * 1.0e18  # um^3 / s
    # C_1 = C * sqrt(16 * gamma_dot * eta / (k_B * T)) in um^(-3/2)
    C1 = C * math.sqrt(16.0 * gamma_dot * eta / (K_B * T)) * 1.0e-9

    def F(a):
        u = C1 * (a ** 1.5)
        if u < 1.0e-5:
            return (a ** 3 / 3.0) * (1.0 - (2.0 / 3.0) * u + 0.5 * u * u)
        return (2.0 / (3.0 * C1 * C1)) * (u - math.log1p(u))

    C0 = F(a0)
    if m >= 0.0:
        target = C0 + m * t
        hi = max(a0 * 2.0, (3.0 * target) ** (1.0 / 3.0) * 2.0, (1.5 * C1 * target) ** (2.0 / 3.0) * 2.0)
        while F(hi) < target:
            hi *= 2.0
        a = brentq(lambda x: F(x) - target, a0, hi, xtol=tol, rtol=tol)
        u = C1 * (a ** 1.5)
        da_dt = (m / (a * a)) * (1.0 + u)
        alpha_eff = (t / a) * da_dt
    else:
        # Unphysical negative m: Eq. (37) shrinks a(t) to 0 in finite time t_collapse = -C0 / m
        t_collapse = -C0 / m
        t_eval = min(t, 0.5 * t_collapse)
        target = C0 + m * t_eval
        a = brentq(lambda x: F(x) - target, 1.0e-12 * a0, a0, xtol=tol, rtol=tol)
        u = C1 * (a ** 1.5)
        da_dt = (m / (a * a)) * (1.0 + u)
        alpha_eff = (t_eval / a) * da_dt

    return {
        "m": float(m),
        "a": float(a),
        "da_dt": float(da_dt),
        "C1_a_three_halves": float(u),
        "alpha_eff": float(alpha_eff),
    }
