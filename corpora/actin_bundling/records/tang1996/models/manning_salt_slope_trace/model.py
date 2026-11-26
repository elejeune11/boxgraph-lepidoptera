def run(inputs, settings):
    return {"alpha": float(inputs["N"])}


def symbolic():
    import sympy as sp

    N = sp.Symbol("N", positive=True)
    return {"alpha": N}
