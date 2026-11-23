def run(inputs, settings):
    return {"n_f": inputs["n_f_0"] * (inputs["D_ratio"] ** 2)}


def symbolic():
    import sympy as sp

    n_f_0, D_ratio = sp.symbols("n_f_0 D_ratio", positive=True)
    return {"n_f": n_f_0 * D_ratio**2}
