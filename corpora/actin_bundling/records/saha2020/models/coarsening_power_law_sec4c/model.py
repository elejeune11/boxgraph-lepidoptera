def run(inputs, settings):
    t = float(inputs["t"])
    beta = float(inputs["beta"])
    return {"R": t ** beta}


def symbolic():
    import sympy as sp

    t, beta = sp.symbols("t beta", positive=True)
    return {"R": t**beta}
