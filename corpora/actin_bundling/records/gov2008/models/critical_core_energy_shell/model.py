def run(inputs, settings):
    return {"E_c_star": inputs["eps_b"] / (inputs["rho"] * inputs["a"] * inputs["lam"])}


def symbolic():
    import sympy as sp

    eps_b, rho, a, lam = sp.symbols("eps_b rho a lam", positive=True)
    return {"E_c_star": eps_b / (rho * a * lam)}
