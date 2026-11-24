def run(inputs, settings):
    return {"m0": (inputs["k_on"] / inputs["k_off"]) * inputs["n0"]}


def symbolic():
    import sympy as sp

    k_off, k_on, n0 = sp.symbols("k_off k_on n0", positive=True)
    return {"m0": (k_on / k_off) * n0}
