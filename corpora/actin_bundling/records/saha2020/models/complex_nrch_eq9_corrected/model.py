def run(inputs, settings):
    c = float(inputs["c"])
    alpha = float(inputs["alpha"])
    kappa = float(inputs["kappa"])
    q = float(inputs["q"])
    q2 = q * q
    return {
        "re_lambda": 4.0 * c * c * q2 - kappa * q2 * q2,
        "im_lambda": alpha * q2,
    }


def symbolic():
    import sympy as sp

    c, alpha, kappa, q = sp.symbols("c alpha kappa q", real=True)
    return {
        "re_lambda": 4 * c**2 * q**2 - kappa * q**4,
        "im_lambda": alpha * q**2,
    }
