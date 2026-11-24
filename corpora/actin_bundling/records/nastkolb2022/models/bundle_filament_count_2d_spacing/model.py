def run(inputs, settings):
    # W in um -> nm is * 1000.0; d_spacing in nm
    return {"n_f": (inputs["W"] * 1000.0) / inputs["d_spacing"]}


def symbolic():
    import sympy as sp

    W, d_spacing = sp.symbols("W d_spacing", positive=True)
    return {"n_f": 1000 * W / d_spacing}
