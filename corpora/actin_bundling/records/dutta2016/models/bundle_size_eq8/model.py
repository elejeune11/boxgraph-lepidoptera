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
    lR = float(inputs["lambda_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    n_k = int(settings.get("n_k", 2000))
    k_max = float(settings.get("k_max", 25.0))
    xtol = float(settings.get("xtol", 1e-10))

    # Check k = 0 macroscopic phase separation (infinite bundle)
    bv0 = 2.0 * math.pi * GR * (lR * lR) - math.pi * GA * (lA * lA)
    if bv0 < 0.0 and rho_star * (-bv0) >= 1.0:
        return {"bundle_size": None, "bundle_size_over_lambda_R": None}

    def g_scalar(k):
        bv = 2.0 * math.pi * GR / (k * k + 1.0 / (lR * lR)) - math.pi * GA * (lA * lA) * math.exp(-0.25 * k * k * lA * lA)
        return rho_star * (-bv) * _b_scalar(k) - 1.0

    ks = np.linspace(0.0, k_max, n_k + 1)
    bv = 2.0 * math.pi * GR / (ks**2 + 1.0 / (lR * lR)) - math.pi * GA * (lA * lA) * np.exp(-0.25 * ks**2 * lA * lA)
    g = rho_star * (-bv) * _b_vec(ks) - 1.0

    pos = np.where(g >= 0.0)[0]
    if len(pos) > 0:
        i = int(pos[0])
        if i == 0:
            return {"bundle_size": None, "bundle_size_over_lambda_R": None}
        k0 = float(brentq(g_scalar, float(ks[i - 1]), float(ks[i]), xtol=xtol))
    else:
        idx = int(np.argmax(g[1:])) + 1
        lo_k = float(ks[max(0, idx - 1)])
        hi_k = float(ks[min(n_k, idx + 1)])
        res = minimize_scalar(lambda x: -g_scalar(x), bounds=(lo_k, hi_k), method="bounded", options={"xatol": xtol})
        if -res.fun < 0.0:
            return {"bundle_size": None, "bundle_size_over_lambda_R": None}
        k0 = float(brentq(g_scalar, lo_k, float(res.x), xtol=xtol))

    if k0 <= 0.0:
        return {"bundle_size": None, "bundle_size_over_lambda_R": None}
    size = 2.0 * math.pi / k0
    return {"bundle_size": size, "bundle_size_over_lambda_R": size / lR}
