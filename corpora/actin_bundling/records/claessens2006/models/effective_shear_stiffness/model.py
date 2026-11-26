def run(inputs, settings):
    return {"k_eff": inputs["k_par"] / inputs["delta"]}


def symbolic():
    import sympy as sp

    k_par, delta = sp.symbols("k_par delta", positive=True)
    return {"k_eff": k_par / delta}
