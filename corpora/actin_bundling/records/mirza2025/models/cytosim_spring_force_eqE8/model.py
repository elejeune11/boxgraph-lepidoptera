def run(inputs, settings):
    return {"f": inputs["K"] * (inputs["u"] - inputs["ell"])}

def symbolic():
    import sympy as sp
    K, u, ell = sp.symbols("K u ell", real=True)
    return {"f": K * (u - ell)}
