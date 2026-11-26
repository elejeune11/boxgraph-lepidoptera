def run(inputs, settings):
    return {"F_actin_poly": inputs["k_B_T"] / inputs["l_mon"]}

def symbolic():
    import sympy as sp
    k_B_T, l_mon = sp.symbols("k_B_T l_mon", positive=True)
    return {"F_actin_poly": k_B_T / l_mon}
