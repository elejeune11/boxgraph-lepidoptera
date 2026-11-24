import math


def run(inputs, settings):
    q = inputs["q"]
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]

    q2 = q * q
    mean_term = -(q2 * (1.0 + kappa) + 2.0 * a + a_Delta + 3.0 * (p1 * p1 + p2 * p2))
    diff_term = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    s = math.sqrt(diff_term * diff_term + 4.0 * rho * rho)
    return {
        "lam_plus": 0.5 * q2 * (mean_term + s),
        "lam_minus": 0.5 * q2 * (mean_term - s),
    }


def symbolic():
    import sympy as sp

    q, a, a_Delta, kappa, phi1_bar, phi2_bar, rho = sp.symbols(
        "q a a_Delta kappa phi1_bar phi2_bar rho", real=True
    )
    mean_term = -(q**2 * (1 + kappa) + 2 * a + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2))
    diff_term = q**2 * (1 - kappa) + 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    s = sp.sqrt(diff_term**2 + 4 * rho**2)
    return {
        "lam_plus": sp.Rational(1, 2) * q**2 * (mean_term + s),
        "lam_minus": sp.Rational(1, 2) * q**2 * (mean_term - s),
    }
