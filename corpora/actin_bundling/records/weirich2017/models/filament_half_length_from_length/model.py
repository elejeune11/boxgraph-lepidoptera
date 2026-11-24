def run(inputs, settings):
    return {"L": 0.5 * inputs["L_f"]}


def symbolic():
    import sympy as sp

    L_f = sp.Symbol("L_f", positive=True)
    return {"L": L_f / 2}
