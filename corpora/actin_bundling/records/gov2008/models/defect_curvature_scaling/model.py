def run(inputs, settings):
    return {"kappa": inputs["a"] / (inputs["L_p"] ** 2)}


def symbolic():
    import sympy as sp

    a, L_p = sp.symbols("a L_p", positive=True)
    return {"kappa": a / L_p**2}
