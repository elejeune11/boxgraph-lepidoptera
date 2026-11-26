def run(inputs, settings):
    return {"b_Kuhn": 2.0 * inputs["l_T"]}


def symbolic():
    import sympy as sp

    l_T = sp.Symbol("l_T", positive=True)
    return {"b_Kuhn": 2 * l_T}
