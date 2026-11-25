def run(inputs, settings):
    return {"n_total": inputs["n_exp"] * inputs["n_per_exp"]}


def symbolic():
    import sympy as sp

    n_exp, n_per_exp = sp.symbols("n_exp n_per_exp", positive=True)
    return {"n_total": n_exp * n_per_exp}
