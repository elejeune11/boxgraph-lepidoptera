def run(inputs, settings):
    return {"omega_imag": inputs["q"] * inputs["c"]}


def symbolic():
    import sympy as sp

    q, c = sp.symbols("q c", positive=True)
    return {"omega_imag": q * c}
