def run(inputs, settings):
    return {"K_d": inputs["k_off"] / inputs["k_on"]}


def symbolic():
    import sympy as sp

    k_off, k_on = sp.symbols("k_off k_on", positive=True)
    return {"K_d": k_off / k_on}
