import math

def run(inputs, settings):
    return {"k_poly": inputs["k_poly_0"] * math.exp(-inputs["F_ext"] / inputs["F_actin_poly"])}

def symbolic():
    import sympy as sp
    k_poly_0, F_ext, F_actin_poly = sp.symbols("k_poly_0 F_ext F_actin_poly", positive=True)
    return {"k_poly": k_poly_0 * sp.exp(-F_ext / F_actin_poly)}
