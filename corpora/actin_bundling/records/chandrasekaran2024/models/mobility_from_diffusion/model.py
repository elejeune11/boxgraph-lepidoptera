def run(inputs, settings):
    # D in um^2/s; kBT in pN*nm -> 1e-3 pN*um => mu in um/(pN*s)
    return {"mu": inputs["D"] / (inputs["kBT"] * 1e-3)}


def symbolic():
    import sympy as sp

    D, kBT = sp.symbols("D kBT", positive=True)
    return {"mu": 1000 * D / kBT}
