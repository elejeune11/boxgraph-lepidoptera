import math

M_TO_NM3 = 0.602214076


def run(inputs, settings):
    Ra = inputs["R_over_a"]
    a_0 = inputs["a_0"]
    a = inputs["a"]
    w = inputs["w"]
    lam0 = 1.0 / inputs["inv_lambda_0"]
    l_B = inputs["l_B"]
    n1 = inputs["n_1"] * M_TO_NM3
    nq = inputs["n_q"] * M_TO_NM3
    q = inputs["q"]

    alpha_b = math.pi * (a * a - a_0 * a_0)
    kappa = math.sqrt(4.0 * math.pi * l_B * (2.0 * n1 + q * (q + 1.0) * nq))
    M = Ra * Ra
    log_term = math.log(lam0 / (q * nq * alpha_b))
    lam_rod = (2.0 / (l_B * q * M)) * log_term
    lam_tot = (M * lam_rod) / (1.0 + kappa * w)
    return {"lambda_tot": lam_tot, "lambda_rod": lam_rod}


def symbolic():
    import sympy as sp

    R_over_a, a_0, a, w, inv_lambda_0, l_B, n_1, n_q, q = sp.symbols(
        "R_over_a a_0 a w inv_lambda_0 l_B n_1 n_q q", positive=True
    )
    alpha_b = sp.pi * (a**2 - a_0**2)
    kappa = sp.sqrt(4 * sp.pi * l_B * M_TO_NM3 * (2 * n_1 + q * (q + 1) * n_q))
    log_term = sp.log((1 / inv_lambda_0) / (q * n_q * M_TO_NM3 * alpha_b))
    lam_rod = 2 * log_term / (l_B * q * R_over_a**2)
    lam_tot = 2 * log_term / (l_B * q * (1 + kappa * w))
    return {"lambda_tot": lam_tot, "lambda_rod": lam_rod}
