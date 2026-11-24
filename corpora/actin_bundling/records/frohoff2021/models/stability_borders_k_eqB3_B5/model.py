import math


def run(inputs, settings):
    q = inputs["k"] / inputs["ell"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    q2 = q * q
    Delta = alpha * alpha - rho * rho
    a_o = -0.5 * (q2 * (1.0 + kappa) + a_Delta + 3.0 * (p1 * p1 + p2 * p2))
    diff = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    disc = diff * diff - 4.0 * Delta
    if disc >= 0.0:
        half_r = 0.5 * math.sqrt(disc)
        a_p = a_o + half_r
        a_m = a_o - half_r
        a_border = a_p
    else:
        a_p = None
        a_m = None
        a_border = a_o
    return {"a_o": a_o, "a_plus": a_p, "a_minus": a_m, "a_border": a_border}


def symbolic():
    import sympy as sp

    k, ell, a_Delta, kappa, phi1_bar, phi2_bar, rho, alpha = sp.symbols(
        "k ell a_Delta kappa phi1_bar phi2_bar rho alpha", real=True
    )
    q = k / ell
    q2 = q**2
    Delta = alpha**2 - rho**2
    a_o = -sp.Rational(1, 2) * (q2 * (1 + kappa) + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2))
    diff = q2 * (1 - kappa) + 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    disc = diff**2 - 4 * Delta
    return {
        "a_o": a_o,
        "a_plus": a_o + sp.Rational(1, 2) * sp.sqrt(disc),
        "a_minus": a_o - sp.Rational(1, 2) * sp.sqrt(disc),
        "a_border": sp.Piecewise((a_o + sp.Rational(1, 2) * sp.sqrt(disc), disc >= 0), (a_o, True)),
    }
