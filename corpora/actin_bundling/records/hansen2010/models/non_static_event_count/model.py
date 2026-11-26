def run(inputs, settings):
    return {"n_remaining": inputs["n_total"] - inputs["n_static"]}


def symbolic():
    import sympy as sp

    n_total, n_static = sp.symbols("n_total n_static", positive=True)
    return {"n_remaining": n_total - n_static}
