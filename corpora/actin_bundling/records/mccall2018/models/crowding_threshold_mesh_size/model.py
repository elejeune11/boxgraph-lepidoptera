K_B = 1.380649e-23  # J/K (exact SI)


def run(inputs, settings):
    xi_m = ((K_B * inputs["T"]) / inputs["Pi_star"]) ** (1.0 / 3.0)
    return {"xi_star": xi_m * 1.0e9}


def symbolic():
    import sympy as sp

    Pi_star, T = sp.symbols("Pi_star T", positive=True)
    return {"xi_star": (sp.Float(K_B * 1.0e27) * T / Pi_star) ** sp.Rational(1, 3)}
