def run(inputs, settings):
    d = inputs["d"]
    return {"F_boundary": inputs["k_boundary"] * d if d > 0.0 else 0.0}


def symbolic():
    import sympy as sp

    k_boundary, d = sp.symbols("k_boundary d", positive=True)
    return {"F_boundary": k_boundary * d}
