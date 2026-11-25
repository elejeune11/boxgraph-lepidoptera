def run(inputs, settings):
    return {"max_error": float(inputs["Nt"]) ** (-1.0)}


def symbolic():
    import sympy as sp

    Nt = sp.Symbol("Nt", positive=True)
    return {"max_error": Nt ** (-1)}
