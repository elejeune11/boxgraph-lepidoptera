def run(inputs, settings):
    return {"R": inputs["c_f"] / inputs["c_a"]}


def symbolic():
    import sympy as sp

    c_f, c_a = sp.symbols("c_f c_a", positive=True)
    return {"R": c_f / c_a}
