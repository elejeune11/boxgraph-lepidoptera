def run(inputs, settings):
    return {"alpha": inputs["l"] / inputs["L"]}


def symbolic():
    import sympy as sp

    l, L = sp.symbols("l L", positive=True)
    return {"alpha": l / L}
