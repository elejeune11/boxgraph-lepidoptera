def run(inputs, settings):
    # From Eq. (8): R_{0,s} = eps_b / (a * rho * k_s_prime), so d(ln R_{0,s}) / d(ln eps_b) = 1 identically.
    _ = (inputs["eps_b"], inputs["k_s_prime"], inputs["a"], inputs["rho"])
    return {"exp_eps_b": 1.0}


def symbolic():
    import sympy as sp

    eps_b, k_s_prime, a, rho = sp.symbols("eps_b k_s_prime a rho", positive=True)
    R0 = eps_b / (a * rho * k_s_prime)
    return {"exp_eps_b": sp.simplify(eps_b * sp.diff(sp.log(R0), eps_b))}
