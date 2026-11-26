def run(inputs, settings):
    return {"k_fil_bind": inputs["k_NMIIA_bind"] * inputs["N_total"]}


def symbolic():
    import sympy as sp
    k_NMIIA_bind, N_total = sp.symbols("k_NMIIA_bind N_total", positive=True)
    return {"k_fil_bind": k_NMIIA_bind * N_total}
