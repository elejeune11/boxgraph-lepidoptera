def run(inputs, settings):
    u_stretch = 0.5 * inputs["k_a"] * ((inputs["r"] - inputs["l_a"]) ** 2)
    u_bend = 0.5 * (inputs["kappa_B"] / inputs["l_a"]) * (inputs["theta"] ** 2)
    return {
        "U_stretch": u_stretch,
        "U_bend": u_bend,
        "U_f": u_stretch + u_bend,
    }


def symbolic():
    import sympy as sp

    k_a, r, l_a, kappa_B = sp.symbols("k_a r l_a kappa_B", positive=True)
    theta = sp.Symbol("theta", real=True)
    u_stretch = sp.Rational(1, 2) * k_a * (r - l_a) ** 2
    u_bend = sp.Rational(1, 2) * (kappa_B / l_a) * theta ** 2
    return {
        "U_stretch": u_stretch,
        "U_bend": u_bend,
        "U_f": u_stretch + u_bend,
    }
