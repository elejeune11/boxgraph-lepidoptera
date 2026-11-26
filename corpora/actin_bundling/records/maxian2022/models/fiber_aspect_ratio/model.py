def run(inputs, settings):
    return {"epsilon": inputs["a"] / inputs["L"]}


def symbolic():
    import sympy as sp

    a, L = sp.symbols("a L", positive=True)
    return {"epsilon": a / L}
