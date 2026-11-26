def run(inputs, settings):
    alpha_s = inputs["alpha_s"]
    a = inputs["a"]
    rho = inputs["rho"]
    r_over_a = 1.0 / (a * a * rho * alpha_s)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}


def symbolic():
    import sympy as sp

    alpha_s, a, rho = sp.symbols("alpha_s a rho", positive=True)
    r_over_a = 1 / (a**2 * rho * alpha_s)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}
