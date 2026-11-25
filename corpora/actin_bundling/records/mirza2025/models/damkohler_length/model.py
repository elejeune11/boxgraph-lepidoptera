import math

def run(inputs, settings):
    return {"ell_D": math.sqrt(inputs["D"] / inputs["k_d"])}

def symbolic():
    import sympy as sp
    D, k_d = sp.symbols("D k_d", positive=True)
    return {"ell_D": sp.sqrt(D / k_d)}
