def run(inputs, settings):
    return {"J_minus": inputs["chi"] * (inputs["k_depoly_minus"] - inputs["k_poly_minus"] * inputs["A_ss"])}

def symbolic():
    import sympy as sp
    chi, k_depoly_minus, k_poly_minus, A_ss = sp.symbols("chi k_depoly_minus k_poly_minus A_ss", positive=True)
    return {"J_minus": chi * (k_depoly_minus - k_poly_minus * A_ss)}
