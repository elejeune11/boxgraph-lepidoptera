def run(inputs, settings):
    return {"K_d_i": inputs["K_d_1"] / (inputs["beta"] ** (inputs["i"] - 1.0))}


def symbolic():
    import sympy as sp

    K_d_1, beta, i = sp.symbols("K_d_1 beta i", positive=True)
    return {"K_d_i": K_d_1 / (beta ** (i - 1))}
