import math


def run(inputs, settings):
    q = inputs["k_l"] ** 2
    alpha = inputs["alpha"]
    P = inputs["P"]
    nu = inputs["nu"]
    beta = inputs["beta"]
    a = 1.0 + beta
    b = 1.0 - beta
    Pi = P / (1.0 + nu * q)
    tr = -q * ((1.0 + alpha) - 2.0 * Pi)
    det = (q**2) * (alpha - Pi * (alpha * a + b))
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        re_lam = 0.5 * (tr + math.sqrt(disc))
        im_lam = 0.0
    else:
        re_lam = 0.5 * tr
        im_lam = 0.5 * math.sqrt(-disc)
    return {
        "tr_tau_L": tr,
        "det_tau_L": det,
        "disc_tau_L": disc,
        "re_lam_plus": re_lam,
        "im_lam_plus": im_lam,
    }


def symbolic():
    import sympy as sp

    k_l, alpha, P, nu = sp.symbols("k_l alpha P nu", positive=True)
    beta = sp.symbols("beta", real=True)
    q = k_l**2
    Pi = P / (1 + nu * q)
    tr = -q * ((1 + alpha) - 2 * Pi)
    det = q**2 * (alpha - Pi * (alpha * (1 + beta) + (1 - beta)))
    disc = sp.expand(tr**2 - 4 * det)
    return {"tr_tau_L": tr, "det_tau_L": det, "disc_tau_L": disc}
