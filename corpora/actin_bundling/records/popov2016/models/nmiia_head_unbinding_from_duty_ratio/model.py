def run(inputs, settings):
    return {"k_unbind_0": inputs["k_bind"] * (1.0 - inputs["rho"]) / inputs["rho"]}

def symbolic():
    import sympy as sp
    k_bind, rho = sp.symbols("k_bind rho", positive=True)
    return {"k_unbind_0": k_bind * (1 - rho) / rho}
