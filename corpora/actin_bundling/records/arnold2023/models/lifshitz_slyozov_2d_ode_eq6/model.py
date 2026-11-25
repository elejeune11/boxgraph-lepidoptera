import math
from scipy.optimize import brentq


def run(inputs, settings):
    sigma = float(inputs["sigma"])
    L = float(inputs["L"])
    a0 = float(inputs["a_0"])
    t = float(inputs["t"])
    tol = float(settings.get("root_tol", 1.0e-10))

    def F(a):
        return (a ** 3 / 9.0) * (1.0 + 3.0 * math.log(L / a))

    C = F(a0)
    target = C + 0.5 * sigma * t
    max_val = L ** 3 / 9.0
    if target >= max_val:
        a = L * (1.0 - 1.0e-9)
    else:
        a = brentq(lambda x: F(x) - target, a0, L * (1.0 - 1.0e-12), xtol=tol, rtol=tol)

    da_dt = -sigma / (2.0 * a * a * math.log(a / L))
    alpha_eff = (t / a) * da_dt
    # Eq. (4) bulk fluctuation evaluated at geometric midpoint r = sqrt(a * L)
    r_mid = math.sqrt(a * L)
    epsilon_mid = (sigma / (2.0 * a)) * (math.log(r_mid / L) / math.log(a / L))
    return {
        "a": float(a),
        "da_dt": float(da_dt),
        "alpha_eff": float(alpha_eff),
        "epsilon_mid": float(epsilon_mid),
    }
