def run(inputs, settings):
    return {"R_c": inputs["K"] / inputs["gamma"]}


def symbolic():
    import sympy as sp

    K, gamma = sp.symbols("K gamma", positive=True)
    return {"R_c": K / gamma}
