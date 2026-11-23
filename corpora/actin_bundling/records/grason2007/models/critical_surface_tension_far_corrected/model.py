import math


def run(inputs, settings):
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    delta_mu = inputs["delta_mu"]
    rho_0 = inputs["rho_0"]
    L = inputs["L"]
    delta = (rho_0 / L) * delta_mu
    sigma_c = (f_0 ** 1.5) / ((3.0 ** 1.5) * (math.pi ** 1.5) * lam3 * math.sqrt(delta))
    return {"Sigma_c": sigma_c}


def symbolic():
    import sympy as sp
    f_0, lambda_3perp, delta_mu, rho_0, L = sp.symbols(
        "f_0 lambda_3perp delta_mu rho_0 L", positive=True
    )
    delta = (rho_0 / L) * delta_mu
    return {
        "Sigma_c": f_0 ** sp.Rational(3, 2)
        / (3 ** sp.Rational(3, 2) * sp.pi ** sp.Rational(3, 2) * lambda_3perp * sp.sqrt(delta))
    }
