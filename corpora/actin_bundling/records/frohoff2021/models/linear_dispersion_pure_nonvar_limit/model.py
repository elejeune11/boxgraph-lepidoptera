import math


def run(inputs, settings):
    q = inputs["q"]
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    alpha = inputs["alpha"]

    q2 = q * q
    mean_term = -(q2 * (1.0 + kappa) + 2.0 * a + a_Delta + 3.0 * (p1 * p1 + p2 * p2))
    diff_term = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    disc = diff_term * diff_term - 4.0 * alpha * alpha
    if disc >= 0.0:
        return {
            "re_lam_plus": 0.5 * q2 * (mean_term + math.sqrt(disc)),
            "im_lam_plus": 0.0,
        }
    return {
        "re_lam_plus": 0.5 * q2 * mean_term,
        "im_lam_plus": 0.5 * q2 * math.sqrt(-disc),
    }


def symbolic():
    import sympy as sp

    q, a, a_Delta, kappa, phi1_bar, phi2_bar, alpha = sp.symbols(
        "q a a_Delta kappa phi1_bar phi2_bar alpha", real=True
    )
    mean_term = -(q**2 * (1 + kappa) + 2 * a + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2))
    diff_term = q**2 * (1 - kappa) + 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    disc = diff_term**2 - 4 * alpha**2
    return {
        "re_lam_plus": sp.Rational(1, 2) * q**2 * (mean_term + sp.sqrt(disc)),
    }
