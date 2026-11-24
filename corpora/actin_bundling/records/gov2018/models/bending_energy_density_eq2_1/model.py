def run(inputs, settings):
    return {
        "f_bend": 0.5 * inputs["kappa"] * (inputs["H"] - inputs["n"] * inputs["H0"]) ** 2
    }


def symbolic():
    import sympy as sp

    kappa, H, n, H0 = sp.symbols("kappa H n H0", real=True)
    return {"f_bend": sp.Rational(1, 2) * kappa * (H - n * H0) ** 2}
