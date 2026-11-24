def run(inputs, settings):
    return {"omega": inputs["gamma_A"] / inputs["gamma_I"]}


def symbolic():
    import sympy as sp

    gamma_A, gamma_I = sp.symbols("gamma_A gamma_I", positive=True)
    return {"omega": gamma_A / gamma_I}
