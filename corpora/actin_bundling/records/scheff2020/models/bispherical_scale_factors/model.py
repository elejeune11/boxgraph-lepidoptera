import math


def run(inputs, settings):
    xi = inputs["xi"]
    eta = inputs["eta"]
    sxi = math.sin(xi)
    Z = 1.0 + sxi * math.cos(eta)
    return {
        "Z": Z,
        "h_phi": sxi * math.sin(eta) / Z,
        "h_xi": 1.0 / Z,
        "h_eta": sxi / Z,
    }


def symbolic():
    import sympy as sp

    xi, eta = sp.symbols("xi eta", positive=True)
    Z = 1 + sp.sin(xi) * sp.cos(eta)
    return {
        "Z": Z,
        "h_phi": sp.sin(xi) * sp.sin(eta) / Z,
        "h_xi": 1 / Z,
        "h_eta": sp.sin(xi) / Z,
    }
