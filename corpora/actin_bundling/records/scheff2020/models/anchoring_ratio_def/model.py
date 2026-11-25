def run(inputs, settings):
    return {"omega": inputs["gamma_A"] / inputs["gamma"]}


def symbolic():
    import sympy as sp

    gamma_A, gamma = sp.symbols("gamma_A gamma", positive=True)
    return {"omega": gamma_A / gamma}
