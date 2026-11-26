def run(inputs, settings):
    return {"N": (inputs["L_0"] * 1000.0) / inputs["l_b"]}


def symbolic():
    import sympy as sp

    L_0, l_b = sp.symbols("L_0 l_b", positive=True)
    return {"N": 1000 * L_0 / l_b}
