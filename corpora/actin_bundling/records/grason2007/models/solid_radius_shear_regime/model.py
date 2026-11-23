def run(inputs, settings):
    return {
        "R_0": ((inputs["mu_perp_bar"] / inputs["Sigma"]) ** (1.0 / 3.0))
        * ((inputs["gamma"] / inputs["mu_par_bar"]) ** (4.0 / 3.0))
    }


def symbolic():
    import sympy as sp
    mu_perp_bar, Sigma, gamma, mu_par_bar = sp.symbols(
        "mu_perp_bar Sigma gamma mu_par_bar", positive=True
    )
    return {
        "R_0": (mu_perp_bar / Sigma) ** sp.Rational(1, 3)
        * (gamma / mu_par_bar) ** sp.Rational(4, 3)
    }
