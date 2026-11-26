import math

def run(inputs, settings):
    L = inputs["L"]
    return {
        "Rg_uniform": 0.5 * L,
        "Rg_max": 0.5 * math.sqrt(3.0) * L,
    }

def symbolic():
    import sympy as sp
    L = sp.symbols("L", positive=True)
    return {
        "Rg_uniform": L / 2,
        "Rg_max": sp.sqrt(3) * L / 2,
    }
