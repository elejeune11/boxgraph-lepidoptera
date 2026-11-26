def run(inputs, settings):
    return {"L_p": inputs["kappa_b"] / inputs["kBT"]}


def symbolic():
    import sympy as sp

    kappa_b, kBT = sp.symbols("kappa_b kBT", positive=True)
    return {"L_p": kappa_b / kBT}
