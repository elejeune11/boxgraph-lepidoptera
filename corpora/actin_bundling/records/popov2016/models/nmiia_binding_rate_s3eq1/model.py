def run(inputs, settings):
    return {"k_fil_bind": inputs["alpha"] * inputs["k_bind"] * inputs["N_total"]}

def symbolic():
    import sympy as sp
    alpha, k_bind, N_total = sp.symbols("alpha k_bind N_total", positive=True)
    return {"k_fil_bind": alpha * k_bind * N_total}
