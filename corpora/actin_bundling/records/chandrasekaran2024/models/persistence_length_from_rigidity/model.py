def run(inputs, settings):
    # kappa in pN*um^2; kBT in pN*nm -> 1e-3 pN*um => l_p in um
    return {"l_p": inputs["kappa"] / (inputs["kBT"] * 1e-3)}


def symbolic():
    import sympy as sp

    kappa, kBT = sp.symbols("kappa kBT", positive=True)
    return {"l_p": 1000 * kappa / kBT}
