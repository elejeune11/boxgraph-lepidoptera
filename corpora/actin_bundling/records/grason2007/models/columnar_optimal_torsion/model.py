def run(inputs, settings):
    R = inputs["R"]
    gamma = inputs["gamma"]
    mu_perp_bar = inputs["mu_perp_bar"]
    lam3 = inputs["lambda_3perp"]
    om0 = (
        (gamma / (4.0 * mu_perp_bar)) ** (1.0 / 3.0)
        * (R ** (-2.0 / 3.0))
        * ((R**2 + lam3**2) ** (-1.0 / 3.0))
    )
    return {"Omega_0": om0}


def symbolic():
    import sympy as sp
    R, gamma, mu_perp_bar, lambda_3perp = sp.symbols(
        "R gamma mu_perp_bar lambda_3perp", positive=True
    )
    return {
        "Omega_0": (gamma / (4 * mu_perp_bar)) ** sp.Rational(1, 3)
        * R ** (-sp.Rational(2, 3))
        * (R**2 + lambda_3perp**2) ** (-sp.Rational(1, 3))
    }
