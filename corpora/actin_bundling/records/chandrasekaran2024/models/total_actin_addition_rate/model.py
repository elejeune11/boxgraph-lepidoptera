def run(inputs, settings):
    return {"R_add": inputs["k_grow"] * inputs["N_fil"]}


def symbolic():
    import sympy as sp

    k_grow, N_fil = sp.symbols("k_grow N_fil", positive=True)
    return {"R_add": k_grow * N_fil}
