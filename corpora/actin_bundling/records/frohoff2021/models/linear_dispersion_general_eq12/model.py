import math


def run(inputs, settings):
    q = inputs["q"]
    f1_pp = inputs["f1_pp"]
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    q2 = q * q
    Delta = alpha * alpha - rho * rho
    tr_B = q2 * (1.0 + Q * kappa) + f1_pp + Q * f2_pp
    det_B = Q * (q2 + f1_pp) * (kappa * q2 + f2_pp) + Q * Delta
    disc = tr_B * tr_B - 4.0 * det_B

    if disc >= 0.0:
        s = math.sqrt(disc)
        lt_p = 0.5 * (-tr_B + s)
        lt_m = 0.5 * (-tr_B - s)
        im_p = 0.0
    else:
        lt_p = -0.5 * tr_B
        lt_m = -0.5 * tr_B
        im_p = 0.5 * q2 * math.sqrt(-disc)

    return {
        "re_lam_plus": q2 * lt_p,
        "re_lam_minus": q2 * lt_m,
        "im_lam_plus": im_p,
        "re_lam_tilde_plus": lt_p,
        "tr_B": tr_B,
        "det_B": det_B,
        "discriminant": disc,
    }


def symbolic():
    import sympy as sp

    q, f1_pp, f2_pp, Q, kappa, rho, alpha = sp.symbols(
        "q f1_pp f2_pp Q kappa rho alpha", real=True
    )
    Delta = alpha**2 - rho**2
    tr_B = q**2 * (1 + Q * kappa) + f1_pp + Q * f2_pp
    det_B = Q * (q**2 + f1_pp) * (kappa * q**2 + f2_pp) + Q * Delta
    disc = tr_B**2 - 4 * det_B
    lt_p = sp.Rational(1, 2) * (-tr_B + sp.sqrt(disc))
    lt_m = sp.Rational(1, 2) * (-tr_B - sp.sqrt(disc))
    return {
        "re_lam_plus": q**2 * lt_p,
        "re_lam_minus": q**2 * lt_m,
        "re_lam_tilde_plus": lt_p,
        "tr_B": tr_B,
        "det_B": det_B,
        "discriminant": disc,
    }
