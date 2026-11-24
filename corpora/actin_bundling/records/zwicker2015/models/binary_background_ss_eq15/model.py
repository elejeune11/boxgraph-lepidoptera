def run(inputs, settings):
    return {"phi_0_B": inputs["k_f"] / (inputs["k_f"] + inputs["k_b"])}


def symbolic():
    import sympy as sp

    k_f, k_b = sp.symbols("k_f k_b", positive=True)
    return {"phi_0_B": k_f / (k_f + k_b)}
