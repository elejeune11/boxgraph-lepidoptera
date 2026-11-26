def run(inputs, settings):
    N = inputs["N"]
    return {"m_max": 2.0 * (N - 1.0) * (3.0 * N - 4.0) / 3.0}


def symbolic():
    import sympy as sp
    N = sp.symbols("N", positive=True)
    return {"m_max": 2 * (N - 1) * (3 * N - 4) / 3}
