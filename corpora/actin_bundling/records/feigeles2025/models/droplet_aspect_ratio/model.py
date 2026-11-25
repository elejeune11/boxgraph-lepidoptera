def run(inputs, settings):
    return {"AR": 2.0 * inputs["R_max"] / inputs["l"]}


def symbolic():
    import sympy as sp
    R_max, l = sp.symbols("R_max l", positive=True)
    return {"AR": 2 * R_max / l}
