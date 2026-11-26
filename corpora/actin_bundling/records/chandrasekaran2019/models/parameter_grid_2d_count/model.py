def run(inputs, settings):
    return {"N_pairs": inputs["N_alpha"] * inputs["N_myosin"]}


def symbolic():
    import sympy as sp
    N_alpha, N_myosin = sp.symbols("N_alpha N_myosin", positive=True)
    return {"N_pairs": N_alpha * N_myosin}
