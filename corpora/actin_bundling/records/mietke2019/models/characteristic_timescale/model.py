def run(inputs, settings):
    return {"tau": inputs["eta_b"] / inputs["gamma"]}


def symbolic():
    import sympy as sp

    eta_b, gamma = sp.symbols("eta_b gamma", positive=True)
    return {"tau": eta_b / gamma}
