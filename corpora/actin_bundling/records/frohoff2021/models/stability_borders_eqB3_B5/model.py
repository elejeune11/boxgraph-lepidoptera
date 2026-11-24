import math


def run(inputs, settings):
    q = inputs["q"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    q2 = q * q
    Delta = alpha * alpha - rho * rho
    a_o = -0.5 * (q2 * (1.0 + kappa) + a_Delta + 3.0 * (p1 * p1 + p2 * p2))
    diff_term = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    disc = diff_term * diff_term - 4.0 * Delta
    if disc < 0.0:
        return {"a_o": a_o, "a_plus": None, "a_minus": None}
    s = 0.5 * math.sqrt(disc)
    return {
        "a_o": a_o,
        "a_plus": a_o + s,
        "a_minus": a_o - s,
    }


def symbolic():
    import sympy as sp

    q, a_Delta, kappa, phi1_bar, phi2_bar, rho, alpha = sp.symbols(
        "q a_Delta kappa phi1_bar phi2_bar rho alpha", real=True
    )
    Delta = alpha**2 - rho**2
    a_o = -sp.Rational(1, 2) * (q**2 * (1 + kappa) + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2))
    diff_term = q**2 * (1 - kappa) + 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    s = sp.Rational(1, 2) * sp.sqrt(diff_term**2 - 4 * Delta)
    return {
        "a_o": a_o,
        "a_plus": a_o + s,
        "a_minus": a_o - s,
    }
