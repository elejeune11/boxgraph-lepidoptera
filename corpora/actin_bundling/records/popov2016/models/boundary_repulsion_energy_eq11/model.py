import math

def run(inputs, settings):
    return {"U_boundary": inputs["eps_boundary"] * math.exp(-inputs["d_i"] / inputs["lambda_boundary"])}

def symbolic():
    import sympy as sp
    eps_boundary, d_i, lambda_boundary = sp.symbols("eps_boundary d_i lambda_boundary", positive=True)
    return {"U_boundary": eps_boundary * sp.exp(-d_i / lambda_boundary)}
