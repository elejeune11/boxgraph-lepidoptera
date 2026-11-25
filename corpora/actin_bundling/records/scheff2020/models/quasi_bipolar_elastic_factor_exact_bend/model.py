import math
import numpy as np


def run(inputs, settings):
    x = float(inputs["x"])
    y = float(inputs["y"])
    n_xi = int(settings.get("n_xi", 48))
    n_eta = int(settings.get("n_eta", 24))
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
    # Exact quasi-bipolar bend term: 4*cos^2(xi) + sin^2(xi)*sin^2(eta)
    num = 4.0 * (cxi * cxi) + (sxi * sxi) * (seta * seta)
    f_el = float(4.0 * math.pi * np.sum(W * h_prod * (t * t) * num / t2N))
    return {"f_el": f_el}
