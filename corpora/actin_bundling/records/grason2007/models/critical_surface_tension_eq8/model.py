import math


def run(inputs, settings):
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    sigma_c = f_0 / ((2.0 ** (4.0 / 3.0)) * math.pi * (lam3 ** (1.0 / 3.0)))
    return {"Sigma_c": sigma_c, "R_0": lam3}


def symbolic():
    import sympy as sp
    f_0, lambda_3perp = sp.symbols("f_0 lambda_3perp", positive=True)
    sigma_c = f_0 / (2 ** sp.Rational(4, 3) * sp.pi * lambda_3perp ** sp.Rational(1, 3))
    return {"Sigma_c": sigma_c, "R_0": lambda_3perp}
