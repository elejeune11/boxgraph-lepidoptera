def run(inputs, settings):
    # W in um -> nm is * 1000.0
    return {"d_spacing": (inputs["W"] * 1000.0) / inputs["n_f"]}


def symbolic():
    import sympy as sp

    W, n_f = sp.symbols("W n_f", positive=True)
    return {"d_spacing": 1000 * W / n_f}
