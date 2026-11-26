def run(inputs, settings):
    return {"L": inputs["N"] * inputs["a"] * 1e-3}


def symbolic():
    import sympy as sp

    N, a = sp.symbols("N a", positive=True)
    return {"L": N * a * sp.Rational(1, 1000)}
