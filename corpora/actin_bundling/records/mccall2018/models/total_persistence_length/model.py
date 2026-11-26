def run(inputs, settings):
    return {"l_T": inputs["l_0"] + inputs["l_E"]}


def symbolic():
    import sympy as sp

    l_0, l_E = sp.symbols("l_0 l_E", positive=True)
    return {"l_T": l_0 + l_E}
