def run(inputs, settings):
    return {"theta": 1.0 - 1.0 / (inputs["N"] * inputs["zeta"])}


def symbolic():
    import sympy as sp

    N, zeta = sp.symbols("N zeta", positive=True)
    return {"theta": 1 - 1 / (N * zeta)}
