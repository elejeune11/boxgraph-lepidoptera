def run(inputs, settings):
    return {"K_ratio": inputs["k_on"] / inputs["k_off"]}


def symbolic():
    import sympy as sp

    k_on, k_off = sp.symbols("k_on k_off", positive=True)
    return {"K_ratio": k_on / k_off}
