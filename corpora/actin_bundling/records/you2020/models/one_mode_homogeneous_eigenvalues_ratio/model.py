import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta_over_kappa"]) * k
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
    return {
        "re_lambda_plus": re_p,
        "re_lambda_minus": re_m,
        "im_lambda_plus": im_p,
    }
