import math
from scipy.optimize import brentq, minimize_scalar


def _omega_cyl(q, p):
    R0 = p["R0"]
    kappa = p["kappa"]
    kBT = p["kBT"]
    Lambda = p["Lambda"]
    H0 = p["H0"]
    ns = p["ns"]
    n0 = p["n0"]
    xi = p["xi"]
    A = p["A"]
    qR2 = (q * R0) ** 2
    H_hh = (kappa / (R0 ** 4)) * (1.0 + qR2 * (1.0 - H0 * R0 * n0) + qR2 * qR2)
    H_hn = (H0 * kappa) / (R0 ** 2)
    H_nn = (H0 * H0 * kappa * n0 + kBT * ns) / n0
    m11 = -H_hh / xi
    m12 = -(H_hn - A) / xi
    mob_q2 = (Lambda * n0 / ns) * (q * q)
    m21 = -mob_q2 * H_hn
    m22 = -mob_q2 * H_nn
    tr = m11 + m22
    det = m11 * m22 - m12 * m21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        return 0.5 * (tr + math.sqrt(disc))
    return 0.5 * tr


def run(inputs, settings):
    R0 = inputs["R0"]
    q_scan_hi = 20.0 / R0
    qs = [q_scan_hi * i / 400.0 for i in range(1, 401)]
    ws = [_omega_cyl(q, inputs) for q in qs]
    best_i = max(range(len(ws)), key=lambda i: ws[i])
    if ws[best_i] <= 0.0:
        return {"q_c": None, "q_max": None, "lambda_max": None, "omega_max": float(ws[best_i])}

    lo = qs[max(0, best_i - 1)]
    hi = qs[min(len(qs) - 1, best_i + 1)]
    res = minimize_scalar(lambda q: -_omega_cyl(q, inputs), bounds=(lo, hi), method="bounded", options={"xatol": 1e-12})
    q_max = float(res.x)
    omega_max = float(-res.fun)

    # Find upper stability boundary q_c > q_max where omega_cyl(q_c) = 0
    q_upper = q_max * 1.5
    while _omega_cyl(q_upper, inputs) > 0.0 and q_upper < 1e5 / R0:
        q_upper *= 2.0
    q_c = float(brentq(lambda q: _omega_cyl(q, inputs), q_max, q_upper, xtol=1e-12))
    return {
        "q_c": q_c,
        "q_max": q_max,
        "lambda_max": 2.0 * math.pi / q_max,
        "omega_max": omega_max,
    }
