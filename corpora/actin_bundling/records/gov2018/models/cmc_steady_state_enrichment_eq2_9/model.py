import math


def run(inputs, settings):
    chi = (inputs["kappa"] / inputs["kBT"]) * (inputs["H"] * inputs["H0"] / inputs["ns"])
    return {"enrichment": math.exp(chi), "chi": chi}


def symbolic():
    import sympy as sp

    kappa, kBT, H, H0, ns = sp.symbols("kappa kBT H H0 ns", real=True)
    chi = (kappa / kBT) * (H * H0 / ns)
    return {"enrichment": sp.exp(chi), "chi": chi}
