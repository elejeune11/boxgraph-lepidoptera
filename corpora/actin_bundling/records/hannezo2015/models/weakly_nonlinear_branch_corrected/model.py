import math


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    phi = float(inputs["phi"])
    i = float(inputs["i"])
    max_rho = float(inputs["max_rho"])
    k = i * math.pi
    chi0 = 1.0 + (k / alpha) ** 2 + phi * (1.0 + (alpha / k) ** 2)
    num = (2.0 * k ** 2 - alpha ** 2) * (alpha ** 2 + k ** 2) * (alpha ** 2 * phi + k ** 2) ** 2
    den = 3.0 * alpha ** 2 * k ** 2 * (4.0 * k ** 4 - alpha ** 4 * phi)
    chi2 = num / den
    eps = max_rho - 1.0
    return {"chi": chi0 + 0.5 * chi2 * eps * eps}


def symbolic():
    import sympy as sp

    alpha, phi, i, max_rho = sp.symbols("alpha phi i max_rho", positive=True)
    k = i * sp.pi
    chi0 = 1 + (k / alpha) ** 2 + phi * (1 + (alpha / k) ** 2)
    chi2 = (2 * k ** 2 - alpha ** 2) * (alpha ** 2 + k ** 2) * (alpha ** 2 * phi + k ** 2) ** 2 / (
        3 * alpha ** 2 * k ** 2 * (4 * k ** 4 - alpha ** 4 * phi)
    )
    return {"chi": chi0 + sp.Rational(1, 2) * chi2 * (max_rho - 1) ** 2}
