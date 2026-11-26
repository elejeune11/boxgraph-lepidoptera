def run(inputs, settings):
    return {"K_minifil": inputs["N_total"] * inputs["K_head"]}

def symbolic():
    import sympy as sp
    N_total, K_head = sp.symbols("N_total K_head", positive=True)
    return {"K_minifil": N_total * K_head}
