def run(inputs, settings):
    return {"ratio": inputs["L_c"] / inputs["epsilon"]}


def symbolic():
    import sympy as sp
    epsilon, L_c = sp.symbols("epsilon L_c", positive=True)
    return {"ratio": L_c / epsilon}
