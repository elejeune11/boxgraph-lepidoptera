def run(inputs, settings):
    return {"D_species": inputs["D_actin"] / inputs["factor"]}

def symbolic():
    import sympy as sp
    D_actin, factor = sp.symbols("D_actin factor", positive=True)
    return {"D_species": D_actin / factor}
