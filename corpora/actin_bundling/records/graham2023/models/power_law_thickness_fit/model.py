def run(inputs, settings):
    return {"T_ring": inputs["a"] * (inputs["D"] ** inputs["b"])}


def symbolic():
    import sympy as sp

    D, a, b = sp.symbols("D a b", positive=True)
    return {"T_ring": a * D ** b}
