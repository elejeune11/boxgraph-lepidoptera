def run(inputs, settings):
    return {"M": inputs["p"] / inputs["s"]}


def symbolic():
    import sympy as sp

    p, s = sp.symbols("p s", positive=True)
    return {"M": p / s}
