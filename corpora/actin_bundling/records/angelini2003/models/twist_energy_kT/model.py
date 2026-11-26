K_B = 1.380649e-23  # J/K (exact SI definition)


def run(inputs, settings):
    d_m = inputs["d"] * 1e-10  # angstrom -> m
    e_j = 0.5 * inputs["kappa"] * (inputs["Delta_theta"] ** 2) / d_m
    return {"E_twist_kT": e_j / (K_B * inputs["T"])}


def symbolic():
    import sympy as sp

    kappa, Delta_theta, d, T = sp.symbols("kappa Delta_theta d T", positive=True)
    return {
        "E_twist_kT": (sp.Rational(1, 2) * kappa * Delta_theta**2 / (d * sp.Rational(1, 10**10)))
        / (sp.Float(K_B, 16) * T)
    }
