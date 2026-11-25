import math
import numpy as np
from scipy.optimize import brentq


def _eval_fn(p1, p2, c11, c12, c21, c22, chi, chip, alpha, kind):
    d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
    d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
    d12 = -((chi + alpha) + 4.0 * p1 * p2 * chip)
    d21 = -((chi - alpha) + 4.0 * p1 * p2 * chip)
    disc = (d11 - d22) ** 2 + 4.0 * d12 * d21
    if kind == "ci":
        return disc
    # C_r boundary: Re(lambda_max) / q^2 = 0
    tr = d11 + d22
    if disc >= 0.0:
        return 0.5 * (tr + math.sqrt(disc))
    return 0.5 * tr


def run(inputs, settings):
    p1 = float(inputs["phi1_bar"])
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])

    kind = str(settings.get("kind", "cr"))
    branch = str(settings.get("branch", "upper"))
    y_lo = float(settings.get("y_lo", -0.6))
    y_hi = float(settings.get("y_hi", 0.6))
    n_scan = int(settings.get("n_scan", 240))

    ys = np.linspace(y_lo, y_hi, n_scan)
    vals = [_eval_fn(p1, float(y), c11, c12, c21, c22, chi, chip, alpha, kind) for y in ys]
    roots = []
    for i in range(len(ys) - 1):
        if vals[i] == 0.0:
            roots.append(float(ys[i]))
        elif vals[i] * vals[i + 1] < 0.0:
            r = brentq(lambda y: _eval_fn(p1, y, c11, c12, c21, c22, chi, chip, alpha, kind), float(ys[i]), float(ys[i + 1]))
            if not roots or abs(r - roots[-1]) > 1e-6:
                roots.append(float(r))
    if not roots:
        return {"phi2_bar": float("nan")}
    if branch == "upper":
        return {"phi2_bar": max(roots)}
    elif branch == "lower":
        return {"phi2_bar": min(roots)}
    elif branch == "second_from_bottom":
        return {"phi2_bar": roots[1] if len(roots) >= 2 else roots[0]}
    return {"phi2_bar": roots[-1]}
