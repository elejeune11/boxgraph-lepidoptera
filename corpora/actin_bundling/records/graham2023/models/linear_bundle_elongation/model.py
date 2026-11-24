def run(inputs, settings):
    return {"y": inputs["v"] * inputs["t"]}


def symbolic():
    import sympy as sp

    v, t = sp.symbols("v t", positive=True)
    return {"y": v * t}
