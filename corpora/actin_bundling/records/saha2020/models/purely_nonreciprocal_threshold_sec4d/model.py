def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    return {"alpha_min": 2.0 * abs(c1 * c1 - c2 * c2)}


def symbolic():
    import sympy as sp

    c1, c2 = sp.symbols("c1 c2", real=True)
    return {"alpha_min": 2 * sp.Abs(c1**2 - c2**2)}
