import math
from scipy.integrate import quad


def run(inputs, settings):
    x = float(inputs["x"])
    Kc = float(inputs["Kc"])
    ell = float(inputs["ell"])
    kBT = float(inputs["kBT"])

    a = math.sqrt(Kc / (2.0 * kBT))
    # Exact analytical integral: int_0^inf x^2 exp(-a^2 (x - ell)^2) dx
    norm = (ell / (2.0 * a * a)) * math.exp(-(a * ell) ** 2) + (
        math.sqrt(math.pi) / (4.0 * a**3)
    ) * (1.0 + 2.0 * (a * ell) ** 2) * math.erfc(-a * ell)
    Z = 1.0 / norm

    def pdf_fn(u):
        return Z * (u * u) * math.exp(-Kc * (u - ell) ** 2 / (2.0 * kBT))

    pdf_val = pdf_fn(x) if x >= 0.0 else 0.0
    x_mode = 0.5 * (ell + math.sqrt(ell * ell + 8.0 * kBT / Kc))

    upper = ell + 12.0 * math.sqrt(kBT / Kc)
    mean_x, _ = quad(lambda u: u * pdf_fn(u), 0.0, upper, epsabs=1e-13, epsrel=1e-12)
    var_x, _ = quad(lambda u: ((u - mean_x) ** 2) * pdf_fn(u), 0.0, upper, epsabs=1e-13, epsrel=1e-12)

    return {
        "pdf": float(pdf_val),
        "x_mode": float(x_mode),
        "mean_x": float(mean_x),
        "std_x": float(math.sqrt(max(0.0, var_x))),
    }


def symbolic():
    import sympy as sp

    x, Kc, ell, kBT = sp.symbols("x Kc ell kBT", positive=True)
    a = sp.sqrt(Kc / (2 * kBT))
    norm = (ell / (2 * a**2)) * sp.exp(-(a * ell) ** 2) + (
        sp.sqrt(sp.pi) / (4 * a**3)
    ) * (1 + 2 * (a * ell) ** 2) * sp.erfc(-a * ell)
    pdf = (x**2 / norm) * sp.exp(-Kc * (x - ell) ** 2 / (2 * kBT))
    x_mode = (ell + sp.sqrt(ell**2 + 8 * kBT / Kc)) / 2
    return {"pdf": pdf, "x_mode": x_mode}
