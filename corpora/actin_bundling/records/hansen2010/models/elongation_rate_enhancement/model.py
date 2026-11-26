def run(inputs, settings):
    return {"fold": inputs["v_vasp"] / inputs["v_basal"]}


def symbolic():
    import sympy as sp

    v_vasp, v_basal = sp.symbols("v_vasp v_basal", positive=True)
    return {"fold": v_vasp / v_basal}
