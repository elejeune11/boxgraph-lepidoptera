def run(inputs, settings):
    return {"mu_perp_bar": (25.0 / 216.0) * inputs["mu_perp"]}


def symbolic():
    import sympy as sp
    mu_perp = sp.Symbol("mu_perp", positive=True)
    return {"mu_perp_bar": sp.Rational(25, 216) * mu_perp}
