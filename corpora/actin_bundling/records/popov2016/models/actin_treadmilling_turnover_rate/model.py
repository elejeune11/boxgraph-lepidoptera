def run(inputs, settings):
    return {"J_tm": inputs["chi"] * (inputs["k_poly_plus"] * inputs["A_c"] - inputs["k_depoly_plus"])}

def symbolic():
    import sympy as sp
    chi, k_poly_plus, k_depoly_plus, A_c = sp.symbols("chi k_poly_plus k_depoly_plus A_c", positive=True)
    return {"J_tm": chi * (k_poly_plus * A_c - k_depoly_plus)}
