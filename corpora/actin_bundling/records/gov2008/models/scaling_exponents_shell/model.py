def run(inputs, settings):
    eps_b = inputs["eps_b"]
    rho = inputs["rho"]
    a = inputs["a"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    ec_star = eps_b / (rho * a * lam)
    if E_c <= ec_star:
        return {"exp_eps_b": None, "exp_rho": None}
    ratio = ec_star / E_c
    val = ratio / (1.0 - ratio)
    return {"exp_eps_b": val, "exp_rho": -val}


def symbolic():
    import sympy as sp

    eps_b, rho, a, lam, E_c = sp.symbols("eps_b rho a lam E_c", positive=True)
    ratio = eps_b / (rho * a * lam * E_c)
    val = ratio / (1 - ratio)
    return {"exp_eps_b": val, "exp_rho": -val}
