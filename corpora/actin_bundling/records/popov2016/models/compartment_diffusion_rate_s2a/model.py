def run(inputs, settings):
    return {"k_diff": inputs["D"] / (inputs["l_comp"] ** 2)}

def symbolic():
    import sympy as sp
    D, l_comp = sp.symbols("D l_comp", positive=True)
    return {"k_diff": D / l_comp**2}
