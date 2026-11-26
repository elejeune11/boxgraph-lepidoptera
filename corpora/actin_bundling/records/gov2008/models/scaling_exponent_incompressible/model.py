def run(inputs, settings):
    # From Eq. (9): R'_{0,b} = a * (pi * (L_p / a)^4 * eps_b / (K * rho))^(1/3), so d(ln R'_{0,b}) / d(ln eps_b) = 1/3 identically.
    _ = (inputs["eps_b"], inputs["K"], inputs["a"], inputs["L_p"], inputs["rho"])
    return {"exp_eps_b": 1.0 / 3.0}


def symbolic():
    import sympy as sp

    eps_b, K, a, L_p, rho = sp.symbols("eps_b K a L_p rho", positive=True)
    R0 = a * (sp.pi * (L_p / a) ** 4 * eps_b / (K * rho)) ** sp.Rational(1, 3)
    return {"exp_eps_b": sp.simplify(eps_b * sp.diff(sp.log(R0), eps_b))}
