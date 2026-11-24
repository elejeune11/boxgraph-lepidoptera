import math


def run(inputs, settings):
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]

    X = 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    Delta = alpha * alpha - rho * rho
    alpha_k0 = math.sqrt(rho * rho + 0.25 * X * X)

    if Delta < 0.0 or abs(1.0 - kappa) < 1e-15:
        return {
            "q_plus_o": None,
            "q_minus_o": None,
            "q_band_lo": None,
            "q_band_hi": None,
            "k_band_lo": None,
            "k_band_hi": None,
            "alpha_k0": alpha_k0,
        }

    two_sd = 2.0 * math.sqrt(Delta)
    qp2 = (two_sd - X) / (1.0 - kappa)
    qm2 = (-two_sd - X) / (1.0 - kappa)
    q_plus_o = math.sqrt(qp2) if qp2 >= 0.0 else None
    q_minus_o = math.sqrt(qm2) if qm2 >= 0.0 else None

    q_lo = q_hi = None
    if X > two_sd and kappa > 1.0:
        # Corrected row 1: since 1 - kappa < 0, 0 < q_+^o < q_-^o
        q_lo, q_hi = q_plus_o, q_minus_o
    elif abs(X) < two_sd:
        if kappa > 1.0:
            q_lo, q_hi = 0.0, q_minus_o
        else:
            q_lo, q_hi = 0.0, q_plus_o
    elif X < -two_sd and kappa < 1.0:
        q_lo, q_hi = q_minus_o, q_plus_o
    elif abs(X - two_sd) <= 1e-15 and kappa > 1.0:
        q_lo, q_hi = q_plus_o, q_minus_o
    elif abs(X + two_sd) <= 1e-15 and kappa < 1.0:
        q_lo, q_hi = q_minus_o, q_plus_o

    return {
        "q_plus_o": q_plus_o,
        "q_minus_o": q_minus_o,
        "q_band_lo": q_lo,
        "q_band_hi": q_hi,
        "k_band_lo": ell * q_lo if q_lo is not None else None,
        "k_band_hi": ell * q_hi if q_hi is not None else None,
        "alpha_k0": alpha_k0,
    }
