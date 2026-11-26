def run(inputs, settings):
    return {"k_on": inputs["K_off"] / inputs["K_d"]}


def symbolic():
    import sympy as sp

    K_off, K_d = sp.symbols("K_off K_d", positive=True)
    return {"k_on": K_off / K_d}
