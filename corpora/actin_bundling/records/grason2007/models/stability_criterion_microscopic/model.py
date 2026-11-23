def run(inputs, settings):
    d = inputs["d"]
    P = inputs["P"]
    lam3 = inputs["lambda_3perp"]
    dE_ratio = inputs["dE_ratio"]
    ratio = ((d / P) ** (4.0 / 3.0)) * ((d / lam3) ** (1.0 / 3.0)) * dE_ratio
    return {"stability_ratio": ratio}


def symbolic():
    import sympy as sp
    d, P, lambda_3perp, dE_ratio = sp.symbols("d P lambda_3perp dE_ratio", positive=True)
    return {
        "stability_ratio": (d / P) ** sp.Rational(4, 3)
        * (d / lambda_3perp) ** sp.Rational(1, 3)
        * dE_ratio
    }
