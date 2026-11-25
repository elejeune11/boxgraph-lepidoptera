import math

def run(inputs, settings):
    return {"ell_s": math.sqrt(inputs["eta"] / inputs["gamma"]) * 1e6}

def symbolic():
    import sympy as sp
    eta, gamma = sp.symbols("eta gamma", positive=True)
    return {"ell_s": sp.sqrt(eta / gamma) * 10**6}
