K_B = 1.380649e-23  # J/K (exact SI definition)


def run(inputs, settings):
    return {"kBT": K_B * inputs["T"]}


def symbolic():
    import sympy as sp

    T = sp.Symbol("T", positive=True)
    return {"kBT": sp.Float(K_B, 16) * T}
