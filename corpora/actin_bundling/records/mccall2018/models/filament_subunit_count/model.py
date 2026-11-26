def run(inputs, settings):
    return {"N_subunits": inputs["L"] / inputs["delta_L"]}


def symbolic():
    import sympy as sp

    L, delta_L = sp.symbols("L delta_L", positive=True)
    return {"N_subunits": L / delta_L}
