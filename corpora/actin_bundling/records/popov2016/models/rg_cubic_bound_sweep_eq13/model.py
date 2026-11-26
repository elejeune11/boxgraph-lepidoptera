import math

def run(inputs, settings):
    L = inputs["L"]
    return {
        "Rg_max": (math.sqrt(3.0) / 2.0) * L,
        "Rg_uniform": 0.5 * L,
    }

def symbolic():
    import sympy as sp
    L, R_alpha_a = sp.symbols("L R_alpha_a", positive=True)
    return {
        "Rg_max": (sp.sqrt(3) / 2) * L + 0 * R_alpha_a,
        "Rg_uniform": L / 2 + 0 * R_alpha_a,
    }
