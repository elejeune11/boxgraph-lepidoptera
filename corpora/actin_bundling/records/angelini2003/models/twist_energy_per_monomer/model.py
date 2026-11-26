def run(inputs, settings):
    d_m = inputs["d"] * 1e-10  # angstrom -> m
    return {"E_twist": 0.5 * inputs["kappa"] * (inputs["Delta_theta"] ** 2) / d_m}


def symbolic():
    import sympy as sp

    kappa, Delta_theta, d = sp.symbols("kappa Delta_theta d", positive=True)
    return {"E_twist": sp.Rational(1, 2) * kappa * Delta_theta**2 / (d * sp.Rational(1, 10**10))}
