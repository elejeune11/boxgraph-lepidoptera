def run(inputs, settings):
    return {"tau": 1.0 / (inputs["k_poly_plus"] * inputs["A"])}

def symbolic():
    import sympy as sp
    k_poly_plus, A = sp.symbols("k_poly_plus A", positive=True)
    return {"tau": 1 / (k_poly_plus * A)}
