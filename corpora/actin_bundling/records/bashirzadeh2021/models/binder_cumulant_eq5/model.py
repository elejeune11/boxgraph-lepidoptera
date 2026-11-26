def run(inputs, settings):
    m2 = inputs["m2"]
    m4 = inputs["m4"]
    return {"U_4": 1.0 - m4 / (3.0 * (m2 ** 2))}


def symbolic():
    import sympy as sp

    m2, m4 = sp.symbols("m2 m4", positive=True)
    return {"U_4": 1 - m4 / (3 * m2 ** 2)}
