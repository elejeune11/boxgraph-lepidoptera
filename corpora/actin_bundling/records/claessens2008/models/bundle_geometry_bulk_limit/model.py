def run(inputs, settings):
    return {"R": 3.0 / inputs["u"]}


def symbolic():
    import sympy as sp

    u = sp.Symbol("u", positive=True)
    return {"R": 3 / u}
