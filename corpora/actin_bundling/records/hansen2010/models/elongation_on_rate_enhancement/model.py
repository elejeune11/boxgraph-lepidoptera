def run(inputs, settings):
    return {"fold": inputs["k_on_vasp"] / inputs["k_on_basal"]}


def symbolic():
    import sympy as sp

    k_on_vasp, k_on_basal = sp.symbols("k_on_vasp k_on_basal", positive=True)
    return {"fold": k_on_vasp / k_on_basal}
