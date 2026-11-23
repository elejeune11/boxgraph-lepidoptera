def run(inputs, settings):
    return {"n_remaining": inputs["n_before"] - inputs["n_excluded"]}


def symbolic():
    import sympy as sp
    n_before, n_excluded = sp.symbols("n_before n_excluded", positive=True)
    return {"n_remaining": n_before - n_excluded}
