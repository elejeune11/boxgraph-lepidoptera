def run(inputs, settings):
    return {"n_af": float(inputs["n_a"]) * float(inputs["p_fa"])}


def symbolic():
    import sympy as sp

    n_a, p_fa = sp.symbols("n_a p_fa", positive=True)
    return {"n_af": n_a * p_fa}
