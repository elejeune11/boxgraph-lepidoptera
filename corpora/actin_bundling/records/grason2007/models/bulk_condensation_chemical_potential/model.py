def run(inputs, settings):
    return {"mu_c": -inputs["epsilon"] * inputs["L"] / inputs["rho_0"]}


def symbolic():
    import sympy as sp
    epsilon, L, rho_0 = sp.symbols("epsilon L rho_0", positive=True)
    return {"mu_c": -epsilon * L / rho_0}
