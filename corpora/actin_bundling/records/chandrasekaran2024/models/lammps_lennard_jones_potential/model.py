def run(inputs, settings):
    eps = inputs["eps"]
    sigma = inputs["sigma"]
    r = inputs["r"]
    r_cut = inputs["r_cut"]
    r_wca = (2.0 ** (1.0 / 6.0)) * sigma
    if r > r_cut:
        u = 0.0
    else:
        sr = sigma / r
        u = 4.0 * eps * (sr**12 - sr**6)
    return {"U": u, "r_wca": r_wca}


def symbolic():
    import sympy as sp

    eps, sigma, r, r_cut = sp.symbols("eps sigma r r_cut", positive=True)
    return {
        "U": 4 * eps * ((sigma / r) ** 12 - (sigma / r) ** 6),
        "r_wca": 2 ** sp.Rational(1, 6) * sigma,
    }
