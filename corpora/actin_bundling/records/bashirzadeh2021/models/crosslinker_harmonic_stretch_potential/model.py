def run(inputs, settings):
    return {"U_xl_stretch": 0.5 * inputs["kappa"] * (inputs["d"] - inputs["l"]) ** 2}


def symbolic():
    import sympy as sp

    kappa, d, l = sp.symbols("kappa d l", real=True)
    return {"U_xl_stretch": sp.Rational(1, 2) * kappa * (d - l) ** 2}
