import math
import numpy as np
from scipy import optimize


def _fs_v(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        fs = math.pi * inv * (8.0 / 3.0 + (8.0 / 5.0) * inv2)
        v = math.pi * inv2 * (16.0 / 15.0 + (32.0 / 105.0) * inv2)
        return fs, v
    x2 = x * x
    x3 = x * x2
    at = math.atan(1.0 / x)
    fs = (2.0 * math.pi / x2) * (x + x3 - (x2 * x2 - 1.0) * at)
    v = (math.pi / (2.0 * x3)) * (x + (2.0 / 3.0) * x3 + x3 * x2 - (x2 - 1.0) * ((x2 + 1.0) ** 2) * at)
    return fs, v


def run(inputs, settings):
    omega = float(inputs["omega"])
    n_xi = int(settings.get("n_xi", 48))
    xtol = float(settings.get("xtol", 1e-8))

    u, wu = np.polynomial.legendre.leggauss(n_xi)
    s = 0.5 * (u + 1.0)
    xi = 0.5 * math.pi * (s * s)
    wxi = math.pi * s * 0.5 * wu
    sxi = np.sin(xi)
    cxi2 = np.cos(xi) ** 2

    def eval_F_hom(x):
        fs, v = _fs_v(x)
        eta0 = 2.0 * math.atan(1.0 / x)
        s0 = math.sin(eta0)
        c0 = math.cos(eta0)
        Z0 = 1.0 + sxi * c0
        # Exact y -> oo limit of Eq. (2C): h_phi*h_xi * cos^2(xi)*sin^2(eta0)/Z0^2
        fan_inf = 4.0 * math.pi * np.sum(wxi * (sxi * (s0**3) * cxi2) / (Z0**4))
        return float((fs + omega * fan_inf) / (v ** (2.0 / 3.0)))

    x_hi = max(10.0, 4.0 * (1.0 + math.sqrt(max(omega, 0.0))))
    res = optimize.minimize_scalar(eval_F_hom, bounds=(1.0, x_hi), method="bounded", options={"xatol": xtol})
    return {"x_eq": float(res.x), "F_nondim": float(res.fun)}
