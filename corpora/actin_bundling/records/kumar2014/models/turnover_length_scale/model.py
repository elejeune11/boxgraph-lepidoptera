import math


def run(inputs, settings):
    ell_k = math.sqrt(inputs["D"] / inputs["kappa"])
    return {"ell_kappa": ell_k, "ell_kappa_over_ell": ell_k / inputs["ell"]}


def symbolic():
    import sympy as sp

    D, kappa, ell = sp.symbols("D kappa ell", positive=True)
    ell_k = sp.sqrt(D / kappa)
    return {"ell_kappa": ell_k, "ell_kappa_over_ell": ell_k / ell}
