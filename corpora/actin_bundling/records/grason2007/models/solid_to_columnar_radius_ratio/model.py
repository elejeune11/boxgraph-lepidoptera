def run(inputs, settings):
    return {"radius_ratio": (inputs["Sigma"] / inputs["gamma"]) ** (8.0 / 3.0)}


def symbolic():
    import sympy as sp
    Sigma, gamma = sp.symbols("Sigma gamma", positive=True)
    return {"radius_ratio": (Sigma / gamma) ** sp.Rational(8, 3)}
