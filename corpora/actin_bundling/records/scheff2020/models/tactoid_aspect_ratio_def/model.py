def run(inputs, settings):
    return {"x": inputs["R"] / inputs["r"]}


def symbolic():
    import sympy as sp

    R, r = sp.symbols("R r", positive=True)
    return {"x": R / r}
