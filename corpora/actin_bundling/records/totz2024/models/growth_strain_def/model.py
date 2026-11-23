def run(inputs, settings):
    return {"g": inputs["dl"] / inputs["l0"]}


def symbolic():
    import sympy as sp
    dl, l0 = sp.symbols("dl l0", positive=True)
    return {"g": dl / l0}
