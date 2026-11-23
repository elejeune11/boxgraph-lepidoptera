import math


def run(inputs, settings):
    # Evaluates the literal algebraic expression printed on p. 3 col 2 top:
    # Sigma_c ~= f_0^(3/2) * 3^(-1/2) * lambda_3perp^(-1) * (mu_c - mu)^(-3/2).
    # Note: this printed formula is dimensionally inconsistent ([f_0^(3/2) * lam3^(-1) * mu^(-3/2)] = L^(-7/2) != E/L^2);
    # see variant critical_surface_tension_far_corrected and questions.yaml.
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    delta_mu = inputs["delta_mu"]
    val = (f_0 ** 1.5) * (3.0 ** -0.5) * (lam3 ** -1.0) * (delta_mu ** -1.5)
    return {"Sigma_c": val}


def symbolic():
    import sympy as sp
    f_0, lambda_3perp, delta_mu = sp.symbols("f_0 lambda_3perp delta_mu", positive=True)
    return {
        "Sigma_c": f_0 ** sp.Rational(3, 2)
        * sp.Integer(3) ** (-sp.Rational(1, 2))
        * lambda_3perp ** (-1)
        * delta_mu ** (-sp.Rational(3, 2))
    }
