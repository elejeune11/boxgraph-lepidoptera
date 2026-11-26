"""Maximum peristaltic mode wavelength: lambda_max = (kappa_f * delta / k_perp)^(1/4)."""


def run(inputs, settings):
    # kappa_f [N*m^2] * delta [m] / k_perp [N/m] = m^4 -> ^(1/4) = m -> * 1e9 = nm
    lam_m = (inputs["kappa_f"] * inputs["delta"] / inputs["k_perp"]) ** 0.25
    return {"lambda_max": lam_m * 1e9}


def symbolic():
    import sympy as sp

    kappa_f, delta, k_perp = sp.symbols("kappa_f delta k_perp", positive=True)
    return {"lambda_max": 1e9 * (kappa_f * delta / k_perp) ** sp.Rational(1, 4)}
