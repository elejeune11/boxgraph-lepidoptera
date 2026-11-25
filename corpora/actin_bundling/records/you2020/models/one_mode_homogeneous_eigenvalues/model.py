import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    mean = -0.5 * (aA + aB)
    disc = (aA - aB) ** 2 + 4.0 * (k * k - d * d)
    if disc >= 0.0:
        half_sq = 0.5 * math.sqrt(disc)
        re_p = mean + half_sq
        re_m = mean - half_sq
        im_p = 0.0
    else:
        re_p = mean
        re_m = mean
        im_p = 0.5 * math.sqrt(-disc)
    delta_ep_over_k = math.sqrt(1.0 + (aA - aB) ** 2 / (4.0 * k * k))
    return {
        "re_lambda_plus": re_p,
        "re_lambda_minus": re_m,
        "im_lambda_plus": im_p,
        "delta_complex_over_kappa": delta_ep_over_k,
    }


def symbolic():
    import sympy as sp

    alpha_A, alpha_B, kappa, delta = sp.symbols("alpha_A alpha_B kappa delta", real=True)
    disc = (alpha_A - alpha_B) ** 2 + 4 * (kappa**2 - delta**2)
    return {
        "re_lambda_plus": -sp.Rational(1, 2) * (alpha_A + alpha_B) + sp.Rational(1, 2) * sp.sqrt(disc),
        "re_lambda_minus": -sp.Rational(1, 2) * (alpha_A + alpha_B) - sp.Rational(1, 2) * sp.sqrt(disc),
        "im_lambda_plus": sp.Rational(1, 2) * sp.sqrt(-disc),
        "delta_complex_over_kappa": sp.sqrt(1 + (alpha_A - alpha_B) ** 2 / (4 * kappa**2)),
    }
