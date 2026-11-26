def run(inputs, settings):
    return {"r": inputs["k_on_F"] / inputs["k_on_alpha"]}


def symbolic():
    import sympy as sp

    k_on_F, k_on_alpha = sp.symbols("k_on_F k_on_alpha", positive=True)
    return {"r": k_on_F / k_on_alpha}
