"""Ridge / actin bundle periodic spacing s = 1 / freq (SI Figs. S10D & S11B)."""


def run(inputs, settings):
    return {"s": 1.0 / float(inputs["freq"])}


def symbolic():
    import sympy as sp

    freq = sp.Symbol("freq", positive=True)
    return {"s": 1 / freq}
