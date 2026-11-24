def run(inputs, settings):
    return {"k_s": inputs["k_f"] + inputs["k"] * inputs["psi"]}


def symbolic():
    import sympy as sp

    k_f, k, psi = sp.symbols("k_f k psi", positive=True)
    return {"k_s": k_f + k * psi}
