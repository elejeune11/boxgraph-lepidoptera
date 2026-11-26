def run(inputs, settings):
    return {"b": 10.0 / inputs["sigma"]}


def symbolic():
    import sympy as sp

    sigma = sp.Symbol("sigma", positive=True)
    return {"b": sp.Integer(10) / sigma}
