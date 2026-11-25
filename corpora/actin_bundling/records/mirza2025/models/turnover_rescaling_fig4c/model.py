def run(inputs, settings):
    return {"k_bar_d": inputs["M"] * inputs["k_bar_d_ref"]}

def symbolic():
    import sympy as sp
    k_bar_d_ref, M = sp.symbols("k_bar_d_ref M", positive=True)
    return {"k_bar_d": M * k_bar_d_ref}
