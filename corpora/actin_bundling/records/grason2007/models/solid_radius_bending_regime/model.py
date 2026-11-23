def run(inputs, settings):
    return {
        "R_0": ((inputs["K3_bar"] / inputs["Sigma"]) ** 0.2)
        * ((inputs["gamma"] / inputs["mu_par_bar"]) ** 0.8)
    }


def symbolic():
    import sympy as sp
    K3_bar, Sigma, gamma, mu_par_bar = sp.symbols(
        "K3_bar Sigma gamma mu_par_bar", positive=True
    )
    return {
        "R_0": (K3_bar / Sigma) ** sp.Rational(1, 5)
        * (gamma / mu_par_bar) ** sp.Rational(4, 5)
    }
