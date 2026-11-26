def run(inputs, settings):
    return {"R_0": inputs["lam"]}


def symbolic():
    import sympy as sp

    lam = sp.symbols("lam", positive=True)
    return {"R_0": lam}
