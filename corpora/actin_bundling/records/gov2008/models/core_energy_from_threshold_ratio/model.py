def run(inputs, settings):
    return {"E_c": inputs["ratio"] * inputs["E_c_star"]}


def symbolic():
    import sympy as sp

    ratio, E_c_star = sp.symbols("ratio E_c_star", positive=True)
    return {"E_c": ratio * E_c_star}
