def run(inputs, settings):
    a = inputs["a"]
    gamma = inputs["gamma"]
    l_B = inputs["l_B"]
    lam_star = inputs["lambda_star_frac"] / inputs["inv_lambda_0"]
    ratio = (7.2 * gamma / (l_B * lam_star * lam_star)) ** (1.0 / 3.0)
    return {"R_eq": a * ratio, "R_eq_over_a": ratio}


def symbolic():
    import sympy as sp

    a, gamma, l_B, inv_lambda_0, lambda_star_frac = sp.symbols(
        "a gamma l_B inv_lambda_0 lambda_star_frac", positive=True
    )
    lam_star = lambda_star_frac / inv_lambda_0
    ratio = (sp.Rational(36, 5) * gamma / (l_B * lam_star**2)) ** sp.Rational(1, 3)
    return {"R_eq": a * ratio, "R_eq_over_a": ratio}
