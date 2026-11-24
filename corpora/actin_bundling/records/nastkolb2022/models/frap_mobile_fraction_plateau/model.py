def run(inputs, settings):
    return {"y": inputs["a"] + inputs["c"]}


def symbolic():
    import sympy as sp

    a, c = sp.symbols("a c", positive=True)
    return {"y": a + c}
