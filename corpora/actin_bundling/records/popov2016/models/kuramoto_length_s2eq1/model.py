import math

def run(inputs, settings):
    return {"l_k": math.sqrt(6.0 * inputs["D"] * inputs["tau"])}

def symbolic():
    import sympy as sp
    D, tau = sp.symbols("D tau", positive=True)
    return {"l_k": sp.sqrt(6 * D * tau)}
