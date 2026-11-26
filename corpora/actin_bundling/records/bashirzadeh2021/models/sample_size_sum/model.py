def run(inputs, settings):
    return {"N_total": inputs["N_a"] + inputs["N_b"]}


def symbolic():
    import sympy as sp

    N_a, N_b = sp.symbols("N_a N_b", positive=True)
    return {"N_total": N_a + N_b}
