def run(inputs, settings):
    a1 = inputs["a1"]
    b1 = inputs["b1"]
    phi = inputs["phi"]
    f = 0.5 * a1 * phi**2 + 0.25 * b1 * phi**4
    mu = a1 * phi + b1 * phi**3
    f_pp = a1 + 3.0 * b1 * phi**2
    return {"f": f, "mu": mu, "f_pp": f_pp}


def symbolic():
    import sympy as sp

    a1, b1, phi = sp.symbols("a1 b1 phi", real=True)
    f = sp.Rational(1, 2) * a1 * phi**2 + sp.Rational(1, 4) * b1 * phi**4
    return {
        "f": f,
        "mu": sp.diff(f, phi),
        "f_pp": sp.diff(f, phi, 2),
    }
