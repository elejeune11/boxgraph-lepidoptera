def run(inputs, settings):
    return {"tau": 1.0 / inputs["k_off"]}


def symbolic():
    import sympy as sp

    k_off = sp.symbols("k_off", positive=True)
    return {"tau": 1 / k_off}
