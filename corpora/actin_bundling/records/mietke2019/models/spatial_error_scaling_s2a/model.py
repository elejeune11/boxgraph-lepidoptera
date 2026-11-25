def run(inputs, settings):
    return {"max_error": float(inputs["N"]) ** (-4.0)}


def symbolic():
    import sympy as sp

    N = sp.Symbol("N", positive=True)
    return {"max_error": N ** (-4)}
