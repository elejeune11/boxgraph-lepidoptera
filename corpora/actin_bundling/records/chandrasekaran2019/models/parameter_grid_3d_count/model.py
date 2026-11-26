def run(inputs, settings):
    return {"N_triads": inputs["N_alpha"] * inputs["N_myosin"] * inputs["N_chi"]}


def symbolic():
    import sympy as sp
    N_alpha, N_myosin, N_chi = sp.symbols("N_alpha N_myosin N_chi", positive=True)
    return {"N_triads": N_alpha * N_myosin * N_chi}
