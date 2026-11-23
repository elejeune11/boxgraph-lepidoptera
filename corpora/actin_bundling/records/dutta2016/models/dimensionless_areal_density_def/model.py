def run(inputs, settings):
    return {"rho_star": inputs["rho"] * (inputs["L"] * inputs["kT"] / inputs["epsilon"])}


def symbolic():
    import sympy as sp
    rho, L, kT, epsilon = sp.symbols("rho L kT epsilon", positive=True)
    return {"rho_star": rho * L * kT / epsilon}
