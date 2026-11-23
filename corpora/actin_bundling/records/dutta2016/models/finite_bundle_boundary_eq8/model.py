import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def _b_vec(k):
    k = np.asarray(k, dtype=float)
    out = np.empty_like(k)
    m = k < 1e-2
    km = k[m]
    out[m] = 1.0 - km**2 / 6.0 + km**4 / 48.0 - km**6 / 480.0
    kn = k[~m]
    out[~m] = (-8.0 + 8.0 * np.exp(-0.5 * kn**2) + 4.0 * kn**2) / (kn**4)
    return out


def _b_scalar(k):
    k = abs(float(k))
    if k < 1e-2:
        k2 = k * k
        return 1.0 - k2 / 6.0 + (k2 * k2) / 48.0 - (k2 * k2 * k2) / 480.0
    return (-8.0 + 8.0 * math.exp(-0.5 * k * k) + 4.0 * k * k) / (k ** 4)


def run(inputs, settings):
    rho_star = float(inputs["rho_star"])
    GR = float(inputs["Gamma_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    n_k = int(settings.get("n_k", 2000))
    k_max = float(settings.get("k_max", 25.0))
    lR_max = float(settings.get("lR_max", 5.0))
    xtol = float(settings.get("xtol", 1e-8))

    num = math.pi * GA * (lA * lA) - 1.0 / rho_star
    if num <= 0.0:
        return {"lambda_R": None}
    lR_inf = math.sqrt(num / (2.0 * math.pi * GR))
    if lR_inf >= lR_max:
        return {"lambda_R": None}

    ks = np.linspace(0.0, k_max, n_k + 1)
    bk = _b_vec(ks)

    def max_g(lR):
        bv = 2.0 * math.pi * GR / (ks**2 + 1.0 / (lR * lR)) - math.pi * GA * (lA * lA) * np.exp(-0.25 * ks**2 * lA * lA)
        g = rho_star * (-bv) * bk - 1.0
        idx = int(np.argmax(g[1:])) + 1
        lo_k = float(ks[max(0, idx - 1)])
        hi_k = float(ks[min(n_k, idx + 1)])

        def neg_g(k):
            bv_s = 2.0 * math.pi * GR / (k * k + 1.0 / (lR * lR)) - math.pi * GA * (lA * lA) * math.exp(-0.25 * k * k * lA * lA)
            return -(rho_star * (-bv_s) * _b_scalar(k) - 1.0)

        res = minimize_scalar(neg_g, bounds=(lo_k, hi_k), method="bounded", options={"xatol": xtol})
        return float(-res.fun)

    lo_lR = lR_inf + 1e-5
    if max_g(lo_lR) <= 0.0:
        return {"lambda_R": None}
    if max_g(lR_max) >= 0.0:
        return {"lambda_R": None}
    lR_sol = float(brentq(max_g, lo_lR, lR_max, xtol=xtol))
    return {"lambda_R": lR_sol}
