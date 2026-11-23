import math


def run(inputs, settings):
    R = inputs["R"]
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    epsilon = inputs["epsilon"]
    Sigma = inputs["Sigma"]
    elast = -f_0 * (R ** (4.0 / 3.0)) / ((lam3**2 + R**2) ** (1.0 / 3.0))
    total = elast - math.pi * (R**2) * epsilon + 2.0 * math.pi * R * Sigma
    return {"dF_over_L": total, "dF_elast_over_L": elast}


def symbolic():
    import sympy as sp
    R, f_0, lambda_3perp, epsilon, Sigma = sp.symbols(
        "R f_0 lambda_3perp epsilon Sigma", positive=True
    )
    elast = -f_0 * R ** sp.Rational(4, 3) / (lambda_3perp**2 + R**2) ** sp.Rational(1, 3)
    total = elast - sp.pi * R**2 * epsilon + 2 * sp.pi * R * Sigma
    return {"dF_over_L": total, "dF_elast_over_L": elast}
