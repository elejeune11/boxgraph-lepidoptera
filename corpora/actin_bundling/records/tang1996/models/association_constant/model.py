def run(inputs, settings):
    return {"K": inputs["c_loc"] / inputs["c_free"]}


def symbolic():
    import sympy as sp

    c_loc, c_free = sp.symbols("c_loc c_free", positive=True)
    return {"K": c_loc / c_free}
