def run(inputs, settings):
    return {"eps_bend": inputs["l_p"] * inputs["k_B_T"] / inputs["l_cyl"]}

def symbolic():
    import sympy as sp
    l_p, k_B_T, l_cyl = sp.symbols("l_p k_B_T l_cyl", positive=True)
    return {"eps_bend": l_p * k_B_T / l_cyl}
