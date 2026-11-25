def run(inputs, settings):
    kappa = float(inputs["kappa"])
    delta = float(inputs["delta"])
    return {"kappa_AB": kappa - delta, "kappa_BA": kappa + delta}


def symbolic():
    import sympy as sp

    kappa, delta = sp.symbols("kappa delta", real=True)
    return {"kappa_AB": kappa - delta, "kappa_BA": kappa + delta}
