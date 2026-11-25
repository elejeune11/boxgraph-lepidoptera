def run(inputs, settings):
    C = float(inputs["C"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])
    return {
        "Pi": Z * C * C - B * C * C * C,
        "Pi_prime": 2.0 * Z * C - 3.0 * B * C * C,
    }


def symbolic():
    import sympy as sp

    C, Z, B = sp.symbols("C Z B", real=True)
    return {
        "Pi": Z * C**2 - B * C**3,
        "Pi_prime": 2 * Z * C - 3 * B * C**2,
    }
