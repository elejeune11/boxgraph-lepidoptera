import math
from scipy.optimize import minimize_scalar


def _omega(q, p):
    q2 = q * q
    q4 = q2 * q2
    D = p["Lambda"] * p["kBT"]
    Gh = (p["kappa"] * q4 + p["sigma"] * q2) / p["xi"]
    Gn = (
        D * q2
        + (p["kappa"] * p["Lambda"] * (p["H0"] ** 2) * p["n0"] / p["ns"]) * q2
        - (2.0 * p["J"] * p["Lambda"] * p["n0"] / p["ns"]) * q2 * (1.0 - q2 / p["ns"])
    )
    m12 = (p["A"] - p["kappa"] * p["H0"] * q2) / p["xi"]
    m21 = -(p["kappa"] * p["Lambda"] * p["H0"] * p["n0"] / p["ns"]) * q4
    tr = -(Gh + Gn)
    det = Gh * Gn - m12 * m21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        return 0.5 * (tr + math.sqrt(disc))
    return 0.5 * tr


def run(inputs, settings):
    # Find upper bracket q_hi where omega(q) < 0
    q_hi = 1.0
    while _omega(q_hi, inputs) > 0.0 and q_hi < 1e5:
        q_hi *= 2.0
    # Scan coarse grid first to confirm positive instability band
    best_q = None
    best_w = 0.0
    n_grid = 200
    for i in range(1, n_grid + 1):
        q = q_hi * i / n_grid
        w = _omega(q, inputs)
        if w > best_w:
            best_w = w
            best_q = q
    if best_q is None:
        return {"q_max": None, "lambda_max": None, "omega_max": 0.0}
    lo = max(1e-9, best_q - q_hi / n_grid)
    hi = best_q + q_hi / n_grid
    res = minimize_scalar(lambda q: -_omega(q, inputs), bounds=(lo, hi), method="bounded", options={"xatol": 1e-12})
    q_max = float(res.x)
    w_max = float(-res.fun)
    return {
        "q_max": q_max,
        "lambda_max": 2.0 * math.pi / q_max,
        "omega_max": w_max,
    }
