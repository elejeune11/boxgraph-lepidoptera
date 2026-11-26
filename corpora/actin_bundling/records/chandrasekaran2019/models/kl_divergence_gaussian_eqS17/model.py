import math


def run(inputs, settings):
    m1, s1 = inputs["mu_1"], inputs["sigma_1"]
    m2, s2 = inputs["mu_2"], inputs["sigma_2"]
    return {
        "D_KL": math.log(s2 / s1) + (s1 * s1 + (m1 - m2) ** 2) / (2.0 * s2 * s2) - 0.5
    }


def symbolic():
    import sympy as sp
    mu_1, sigma_1, mu_2, sigma_2 = sp.symbols("mu_1 sigma_1 mu_2 sigma_2", positive=True)
    return {
        "D_KL": sp.log(sigma_2 / sigma_1)
        + (sigma_1 ** 2 + (mu_1 - mu_2) ** 2) / (2 * sigma_2 ** 2)
        - sp.Rational(1, 2)
    }
