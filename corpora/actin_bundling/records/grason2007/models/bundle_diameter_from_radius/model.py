def run(inputs, settings):
    return {"D": 2.0 * inputs["R"]}


def symbolic():
    import sympy as sp
    R = sp.Symbol("R", positive=True)
    return {"D": 2 * R}
