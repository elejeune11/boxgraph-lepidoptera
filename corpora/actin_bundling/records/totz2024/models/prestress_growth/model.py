def run(inputs, settings):
    return {"sigma_A": inputs["Y"] * inputs["g"]}


def symbolic():
    import sympy as sp
    Y, g = sp.symbols("Y g", positive=True)
    return {"sigma_A": Y * g}
