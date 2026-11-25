import math
import numpy as np


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
    V = float(inputs["V"])
    x = float(inputs["x"])
    y = float(inputs["y"])
    Kg = float(inputs["K_over_gamma"])
    omega = float(inputs["omega"])
    n_xi = int(settings.get("n_xi", 48))
    n_eta = int(settings.get("n_eta", 24))

    fs, v = _fs_v(x)
    eta0 = 2.0 * math.atan(1.0 / x)
    t = 1.0 / (y * y)
    u, wu = np.polynomial.legendre.leggauss(n_xi)
    s = 0.5 * (u + 1.0)
    xi = 0.5 * math.pi * (s * s)
    wxi = math.pi * s * 0.5 * wu
    sxi = np.sin(xi)
    cxi = np.cos(xi)
    sxi2 = sxi * sxi
    cxi2 = cxi * cxi

    vn, wv = np.polynomial.legendre.leggauss(n_eta)
    eta = 0.5 * eta0 * (vn + 1.0)
    weta = 0.5 * eta0 * wv
    seta = np.sin(eta)[None, :]
    ceta = np.cos(eta)[None, :]
    Z = 1.0 + sxi[:, None] * ceta
    t2N = (t * sxi[:, None] * ceta + 0.5 * (1.0 - t) * Z) ** 2 + t * sxi2[:, None] * (seta * seta)
    h_prod = (sxi2[:, None] * seta) / (Z * Z * Z)
    num = 4.0 * cxi2[:, None] + sxi2[:, None] * (ceta * ceta)
    W = wxi[:, None] * weta[None, :]
    fel = float(4.0 * math.pi * np.sum(W * h_prod * (t * t) * num / t2N))

    s0 = math.sin(eta0)
    c0 = math.cos(eta0)
    Z0 = 1.0 + sxi * c0
    t2N0 = (t * sxi * c0 + 0.5 * (1.0 - t) * Z0) ** 2 + t * sxi2 * (s0 * s0)
    fan = float(4.0 * math.pi * np.sum(wxi * (sxi * s0 / (Z0 * Z0)) * (((1.0 - t) ** 2) * cxi2 * (s0 * s0)) / (4.0 * t2N0)))

    c = Kg / (V ** (1.0 / 3.0))
    F_nondim = c * fel / (v ** (1.0 / 3.0)) + (fs + omega * fan) / (v ** (2.0 / 3.0))
    return {"F_nondim": F_nondim}
