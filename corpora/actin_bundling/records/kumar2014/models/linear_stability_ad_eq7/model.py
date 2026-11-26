import math


def run(inputs, settings):
    q = inputs["k_l"] ** 2
    alpha = inputs["alpha"]
    P = inputs["P"]
    nu = inputs["nu"]
    a = inputs["A0_fA"]
    b = inputs["I0_fI"]
    Pi = P / (1.0 + nu * q)
    # Matrix entries of tau*L:
    # M11 = -q + q*Pi*a, M22 = -alpha*q + q*Pi*b, M12*M21 = (q*Pi)^2 * a * b
    M11 = q * (-1.0 + Pi * a)
    M22 = q * (-alpha + Pi * b)
    tr = M11 + M22
    det = M11 * M22 - (q * Pi) ** 2 * a * b
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
    A0_fA, I0_fI = sp.symbols("A0_fA I0_fI", real=True)
    q = k_l**2
    Pi = P / (1 + nu * q)
    M11 = q * (-1 + Pi * A0_fA)
    M22 = q * (-alpha + Pi * I0_fI)
    tr = M11 + M22
    det = sp.simplify(M11 * M22 - (q * Pi) ** 2 * A0_fA * I0_fI)
    disc = sp.expand(tr**2 - 4 * det)
    return {
        "tr_tau_L": tr,
        "det_tau_L": det,
        "disc_tau_L": disc,
    }
