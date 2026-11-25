def run(inputs, settings):
    return {"phi_c": inputs["D"] / inputs["L"]}


def symbolic():
    import sympy as sp

    L, D = sp.symbols("L D", positive=True)
    return {"phi_c": D / L}
