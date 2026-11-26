K_B = 1.380649e-23  # J/K (exact SI Boltzmann constant)


def run(inputs, settings):
    return {"kappa_a": inputs["l_p"] * K_B * inputs["T"]}


def symbolic():
    import sympy as sp

    l_p, T = sp.symbols("l_p T", positive=True)
    return {"kappa_a": l_p * K_B * T}
