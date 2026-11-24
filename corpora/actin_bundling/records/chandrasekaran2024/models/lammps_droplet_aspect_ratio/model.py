def run(inputs, settings):
    return {"aspect_ratio": inputs["c"] / inputs["a"]}


def symbolic():
    import sympy as sp

    a, c = sp.symbols("a c", positive=True)
    return {"aspect_ratio": c / a}
