def run(inputs, settings):
    c_free = inputs["c_free"]
    return {"c_total": c_free + inputs["k"] * inputs["v"] * c_free}


def symbolic():
    import sympy as sp

    c_free, k, v = sp.symbols("c_free k v", positive=True)
    return {"c_total": c_free + k * v * c_free}
