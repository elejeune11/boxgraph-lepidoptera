def run(inputs, settings):
    return {"ratio": inputs["epsilon"] / inputs["L_c"]}


def symbolic():
    import sympy as sp
    epsilon, L_c = sp.symbols("epsilon L_c", positive=True)
    return {"ratio": epsilon / L_c}
