def run(inputs, settings):
    return {"L": inputs["N_mon"] / inputs["rho_mon"]}


def symbolic():
    import sympy as sp

    N_mon, rho_mon = sp.symbols("N_mon rho_mon", positive=True)
    return {"L": N_mon / rho_mon}
