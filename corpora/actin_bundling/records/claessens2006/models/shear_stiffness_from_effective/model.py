def run(inputs, settings):
    return {"k_par": inputs["k_eff"] * inputs["delta"]}


def symbolic():
    import sympy as sp

    k_eff, delta = sp.symbols("k_eff delta", positive=True)
    return {"k_par": k_eff * delta}
