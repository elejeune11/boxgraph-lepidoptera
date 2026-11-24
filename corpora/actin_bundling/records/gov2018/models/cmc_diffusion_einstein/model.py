def run(inputs, settings):
    return {"D": inputs["Lambda"] * inputs["kBT"]}


def symbolic():
    import sympy as sp

    Lambda, kBT = sp.symbols("Lambda kBT", positive=True)
    return {"D": Lambda * kBT}
