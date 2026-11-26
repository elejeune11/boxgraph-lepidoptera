def run(inputs, settings):
    num = inputs["k_depoly_plus"] + inputs["k_depoly_minus"]
    den = inputs["k_poly_plus"] + inputs["k_poly_minus"]
    return {"A_c": num / den}

def symbolic():
    import sympy as sp
    k_depoly_plus, k_depoly_minus, k_poly_plus, k_poly_minus = sp.symbols(
        "k_depoly_plus k_depoly_minus k_poly_plus k_poly_minus", positive=True
    )
    return {"A_c": (k_depoly_plus + k_depoly_minus) / (k_poly_plus + k_poly_minus)}
