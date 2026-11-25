import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    kd = k * k - d * d
    tr = 2.0 * aA - aB - 3.0 * kd / aB
    diff = 2.0 * aA + aB - 3.0 * kd / aB
    disc = diff * diff + 4.0 * kd
    if disc >= 0.0:
        sq = math.sqrt(disc)
        re_p = 0.5 * (tr + sq)
        re_m = 0.5 * (tr - sq)
    else:
        re_p = 0.5 * tr
        re_m = 0.5 * tr
    dc2 = k * k + aB * aB
    lam_theta = -(dc2 - d * d) / aB
    return {
        "re_lambda_s_plus": re_p,
        "re_lambda_s_minus": re_m,
        "lambda_theta_s": lam_theta,
    }


def symbolic():
    import sympy as sp

    alpha_A, alpha_B, kappa, delta = sp.symbols("alpha_A alpha_B kappa delta", real=True)
    kd = kappa**2 - delta**2
    tr = 2 * alpha_A - alpha_B - 3 * kd / alpha_B
    diff = 2 * alpha_A + alpha_B - 3 * kd / alpha_B
    disc = diff**2 + 4 * kd
    dc2 = kappa**2 + alpha_B**2
    return {
        "re_lambda_s_plus": sp.Rational(1, 2) * tr + sp.Rational(1, 2) * sp.sqrt(disc),
        "re_lambda_s_minus": sp.Rational(1, 2) * tr - sp.Rational(1, 2) * sp.sqrt(disc),
        "lambda_theta_s": -(dc2 - delta**2) / alpha_B,
    }
