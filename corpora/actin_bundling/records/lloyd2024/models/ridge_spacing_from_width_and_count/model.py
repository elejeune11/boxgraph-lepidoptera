"""Average longitudinal ridge spacing s = W / N (SI p. 5)."""


def run(inputs, settings):
    return {"s": float(inputs["W"]) / float(inputs["N"])}


def symbolic():
    import sympy as sp

    W, N = sp.symbols("W N", positive=True)
    return {"s": W / N}
