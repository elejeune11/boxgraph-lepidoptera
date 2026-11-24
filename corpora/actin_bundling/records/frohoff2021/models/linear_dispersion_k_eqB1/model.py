import math


def run(inputs, settings):
    q = inputs["k"] / inputs["ell"]
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    q2 = q * q
    Delta = alpha * alpha - rho * rho
    s_sum = q2 * (1.0 + kappa) + 2.0 * a + a_Delta + 3.0 * (p1 * p1 + p2 * p2)
    diff = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    disc = diff * diff - 4.0 * Delta
    if disc >= 0.0:
        r = math.sqrt(disc)
        re_p = 0.5 * q2 * (-s_sum + r)
        re_m = 0.5 * q2 * (-s_sum - r)
        im_p = 0.0
    else:
        re_p = -0.5 * q2 * s_sum
        re_m = re_p
        im_p = 0.5 * q2 * math.sqrt(-disc)
    return {"re_lam_plus": re_p, "re_lam_minus": re_m, "im_lam_plus": im_p}


def symbolic():
    import sympy as sp

    k, ell, a, a_Delta, kappa, phi1_bar, phi2_bar, rho, alpha = sp.symbols(
        "k ell a a_Delta kappa phi1_bar phi2_bar rho alpha", real=True
    )
    q = k / ell
    q2 = q**2
    Delta = alpha**2 - rho**2
    s_sum = q2 * (1 + kappa) + 2 * a + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2)
    diff = q2 * (1 - kappa) + 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    disc = diff**2 - 4 * Delta
    return {
        "re_lam_plus": sp.Piecewise(
            (sp.Rational(1, 2) * q2 * (-s_sum + sp.sqrt(disc)), disc >= 0),
            (-sp.Rational(1, 2) * q2 * s_sum, True),
        ),
        "re_lam_minus": sp.Piecewise(
            (sp.Rational(1, 2) * q2 * (-s_sum - sp.sqrt(disc)), disc >= 0),
            (-sp.Rational(1, 2) * q2 * s_sum, True),
        ),
        "im_lam_plus": sp.Piecewise(
            (sp.Integer(0), disc >= 0),
            (sp.Rational(1, 2) * q2 * sp.sqrt(-disc), True),
        ),
    }
