import math


def run(inputs, settings):
    q2 = inputs["k_l"] ** 2
    d = inputs["Da_over_Ds"]
    rho_a = inputs["rho_a"]
    rho_s = inputs["rho_s"]
    a11 = -d * q2 + rho_a
    a12 = rho_a
    a21 = -2.0 * rho_s
    a22 = -q2 - rho_s
    tr = a11 + a22
    det = a11 * a22 - a12 * a21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        re_lam = 0.5 * (tr + math.sqrt(disc))
        im_lam = 0.0
    else:
        re_lam = 0.5 * tr
        im_lam = 0.5 * math.sqrt(-disc)
    return {
        "re_lam_dimless": float(re_lam),
        "im_lam_dimless": float(im_lam),
        "tr_dimless": float(tr),
        "det_dimless": float(det),
    }


def symbolic():
    import sympy as sp

    k_l, Da_over_Ds, rho_a, rho_s = sp.symbols("k_l Da_over_Ds rho_a rho_s", positive=True)
    q2 = k_l**2
    a11 = -Da_over_Ds * q2 + rho_a
    a12 = rho_a
    a21 = -2 * rho_s
    a22 = -q2 - rho_s
    tr = a11 + a22
    det = a11 * a22 - a12 * a21
    disc = tr**2 - 4 * det
    re_lam = sp.Piecewise((sp.Rational(1, 2) * (tr + sp.sqrt(disc)), disc >= 0), (sp.Rational(1, 2) * tr, True))
    im_lam = sp.Piecewise((sp.Integer(0), disc >= 0), (sp.Rational(1, 2) * sp.sqrt(-disc), True))
    return {
        "re_lam_dimless": re_lam,
        "im_lam_dimless": im_lam,
        "tr_dimless": tr,
        "det_dimless": det,
    }
