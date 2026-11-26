def run(inputs, settings):
    return {"R": inputs["D"] / 2.0}


def symbolic():
    import sympy as sp

    D = sp.Symbol("D", positive=True)
    return {"R": D / 2}
