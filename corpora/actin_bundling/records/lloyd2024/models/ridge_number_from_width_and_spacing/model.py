"""Total ridge number N = W / s, algebraic inverse of SI p. 5 definition s = W / N."""


def run(inputs, settings):
    return {"N": float(inputs["W"]) / float(inputs["s"])}


def symbolic():
    import sympy as sp

    W, s = sp.symbols("W s", positive=True)
    return {"N": W / s}
