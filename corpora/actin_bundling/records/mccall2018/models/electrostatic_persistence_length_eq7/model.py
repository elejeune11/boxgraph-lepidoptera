def run(inputs, settings):
    return {"l_E": (inputs["l_B"] * (inputs["f"] ** 2)) / (4.0 * inputs["kappa_sq"] * (inputs["b"] ** 2))}


def symbolic():
    import sympy as sp

    l_B, f, kappa_sq, b = sp.symbols("l_B f kappa_sq b", positive=True)
    return {"l_E": l_B * f ** 2 / (4 * kappa_sq * b ** 2)}
