def run(inputs, settings):
    return {"K_b_ring": inputs["K_b_fil"] * inputs["N_fil"]}


def symbolic():
    import sympy as sp

    K_b_fil, N_fil = sp.symbols("K_b_fil N_fil", positive=True)
    return {"K_b_ring": K_b_fil * N_fil}
