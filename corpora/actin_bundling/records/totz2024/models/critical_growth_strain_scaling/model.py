def run(inputs, settings):
    return {"g_c": inputs["L"] ** (-2)}


def symbolic():
    import sympy as sp
    L = sp.symbols("L", positive=True)
    return {"g_c": L**(-2)}
