def run(inputs, settings):
    return {"kappa": inputs["K_d"] / inputs["c_a"]}


def symbolic():
    import sympy as sp

    K_d, c_a = sp.symbols("K_d c_a", positive=True)
    return {"kappa": K_d / c_a}
