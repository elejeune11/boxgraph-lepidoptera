def run(inputs, settings):
    return {"l_p": inputs["kappa"] / inputs["kBT"]}


def symbolic():
    import sympy as sp

    kappa, kBT = sp.symbols("kappa kBT", positive=True)
    return {"l_p": kappa / kBT}
