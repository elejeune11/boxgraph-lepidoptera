import math


def run(inputs, settings):
    boltz = math.exp(-0.5 * inputs["Kc"] * (inputs["ell_k"] - inputs["ell"]) ** 2 / inputs["kBT"])
    kons = inputs["kons0"] * boltz
    return {
        "kons": kons,
        "rate_ratio": kons / inputs["koffs0"],
    }


def symbolic():
    import sympy as sp

    kons0, koffs0, Kc, ell, ell_k, kBT = sp.symbols("kons0 koffs0 Kc ell ell_k kBT", positive=True)
    boltz = sp.exp(-Kc * (ell_k - ell) ** 2 / (2 * kBT))
    return {
        "kons": kons0 * boltz,
        "rate_ratio": (kons0 / koffs0) * boltz,
    }
