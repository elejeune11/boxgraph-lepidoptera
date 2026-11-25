import math
import numpy as np


def _fs(x):
    if x > 1e3:
        inv = 1.0 / x
        return math.pi * inv * (8.0 / 3.0 + (8.0 / 5.0) * inv * inv)
    x2 = x * x
    return (2.0 * math.pi / x2) * (x + x * x2 - (x2 * x2 - 1.0) * math.atan(1.0 / x))


def _fel_fan_printed(x, y, n_xi, n_eta):
    eta0 = 2.0 * math.atan(1.0 / x)
    t = 1.0 / (y * y)
    u, wu = np.polynomial.legendre.leggauss(n_xi)
    s = 0.5 * (u + 1.0)
    ws = 0.5 * wu
    xi = 0.5 * math.pi * (s * s)
    wxi = math.pi * s * ws
    v_n, wv = np.polynomial.legendre.leggauss(n_eta)
    eta = 0.5 * eta0 * (v_n + 1.0)
    weta = 0.5 * eta0 * wv

    XI = xi[:, None]
    ETA = eta[None, :]
    W = wxi[:, None] * weta[None, :]
    sxi = np.sin(XI)
    cxi = np.cos(XI)
    seta = np.sin(ETA)
    ceta = np.cos(ETA)
    Z = 1.0 + sxi * ceta
    t2N = (t * sxi * ceta + 0.5 * (1.0 - t) * Z) ** 2 + t * (sxi * sxi) * (seta * seta)
    h_prod = (sxi * sxi * seta) / (Z * Z * Z)
    num = 4.0 * (cxi * cxi) + (sxi * sxi) * (ceta * ceta)
    f_el = float(4.0 * math.pi * np.sum(W * h_prod * (t * t) * num / t2N))

    s0 = math.sin(eta0)
    c0 = math.cos(eta0)
    sxi1 = np.sin(xi)
    cxi1 = np.cos(xi)
    Z0 = 1.0 + sxi1 * c0
    t2N0 = (t * sxi1 * c0 + 0.5 * (1.0 - t) * Z0) ** 2 + t * (sxi1 * sxi1) * (s0 * s0)
    h_surf = (sxi1 * s0) / (Z0 * Z0)
    f_an = float(4.0 * math.pi * np.sum(wxi * h_surf * (((1.0 - t) ** 2) * (cxi1 * cxi1) * (s0 * s0)) / (4.0 * t2N0)))
    return f_el, f_an


def run(inputs, settings):
    R = float(inputs["R"])
    x = float(inputs["x"])
    y = float(inputs["y"])
    K = float(inputs["K"])
    gamma = float(inputs["gamma"])
    omega = float(inputs["omega"])
    n_xi = int(settings.get("n_xi", 48))
    n_eta = int(settings.get("n_eta", 24))
    f_el, f_an = _fel_fan_printed(x, y, n_xi, n_eta)
    F = K * R * f_el + gamma * R * R * (_fs(x) + omega * f_an)
    return {"F": F}
