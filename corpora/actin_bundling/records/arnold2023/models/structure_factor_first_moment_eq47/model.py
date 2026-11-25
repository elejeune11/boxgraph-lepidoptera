"""Compute domain size a(t) and peak wave vector k_max from a static structure factor curve S(k) via SI Eq. (47)."""

import numpy as np


def run(inputs, settings):
    curve = inputs["S_k"]
    ks = np.asarray(curve["x"], dtype=np.float64)
    sk = np.asarray(curve["y"], dtype=np.float64)
    order = np.argsort(ks)
    ks, sk = ks[order], sk[order]

    k_max = float(ks[int(np.argmax(sk))])
    sel = ks <= 2.0 * k_max + 1.0e-12
    num = float(np.trapezoid(sk[sel], ks[sel]))
    den = float(np.trapezoid(ks[sel] * sk[sel], ks[sel]))
    a_sf = num / den
    return {"a_sf": a_sf, "k_max": k_max}
