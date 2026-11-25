def run(inputs, settings):
    return {"t_0": (inputs["ell_s"] ** 2) / inputs["D"]}

def symbolic():
    import sympy as sp
    ell_s, D = sp.symbols("ell_s D", positive=True)
    return {"t_0": ell_s**2 / D}
