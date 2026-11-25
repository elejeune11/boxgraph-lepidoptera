def run(inputs, settings):
    return {"rae": inputs["rho_core"] / inputs["rho_periphery"]}


def symbolic():
    import sympy as sp

    rho_core, rho_periphery = sp.symbols("rho_core rho_periphery", positive=True)
    return {"rae": rho_core / rho_periphery}
