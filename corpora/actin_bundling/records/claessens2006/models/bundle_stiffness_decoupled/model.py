def run(inputs, settings):
    return {"kappa_B": inputs["n"] * inputs["kappa_a"]}


def symbolic():
    import sympy as sp

    n, kappa_a = sp.symbols("n kappa_a", positive=True)
    return {"kappa_B": n * kappa_a}
