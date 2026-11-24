import math


def run(inputs, settings):
    q = inputs["q"]
    L0 = inputs["L0"]
    tau0 = inputs["tau0"]
    f1_pp = inputs["f1_pp"]
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    qF = q * L0
    q2 = qF * qF
    Delta = alpha * alpha - rho * rho
    tr_B = q2 * (1.0 + Q * kappa) + f1_pp + Q * f2_pp
    det_B = Q * (q2 + f1_pp) * (kappa * q2 + f2_pp) + Q * Delta
    disc = tr_B * tr_B - 4.0 * det_B

    if disc >= 0.0:
        s = math.sqrt(disc)
        if tr_B > 0.0:
            lt_p = -2.0 * det_B / (tr_B + s)
            lt_m = -0.5 * (tr_B + s)
        else:
            lt_p = 0.5 * (-tr_B + s)
            denom = -tr_B + s
            lt_m = 2.0 * det_B / denom if denom > 0.0 else 0.0
        im_p = 0.0
    else:
        lt_p = -0.5 * tr_B
        lt_m = -0.5 * tr_B
        im_p = 0.5 * q2 * math.sqrt(-disc)

    return {
        "re_lam_plus": (q2 * lt_p) / tau0,
        "re_lam_minus": (q2 * lt_m) / tau0,
        "im_lam_plus": im_p / tau0,
    }


def symbolic():
    import sympy as sp

    q, L0, tau0, f1_pp, f2_pp, Q, kappa, rho, alpha = sp.symbols(
        "q L0 tau0 f1_pp f2_pp Q kappa rho alpha", real=True
    )
    q2 = (q * L0) ** 2
    Delta = alpha**2 - rho**2
    tr_B = q2 * (1 + Q * kappa) + f1_pp + Q * f2_pp
    det_B = Q * (q2 + f1_pp) * (kappa * q2 + f2_pp) + Q * Delta
    disc = tr_B**2 - 4 * det_B
    lt_p = sp.Rational(1, 2) * (-tr_B + sp.sqrt(disc))
    lt_m = sp.Rational(1, 2) * (-tr_B - sp.sqrt(disc))
    return {
        "re_lam_plus": (q2 * lt_p) / tau0,
        "re_lam_minus": (q2 * lt_m) / tau0,
    }
