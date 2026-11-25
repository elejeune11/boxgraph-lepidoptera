import math


def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])
    alpha = float(inputs["alpha"])
    rad = alpha * alpha - chi * chi
    if rad < 0.0:
        return {"hyperbola_offset": float("nan"), "hyperbola_offset_minus": float("nan")}
    s = math.sqrt(rad) / 6.0
    mid = (c1 * c1 - c2 * c2) / 3.0
    return {
        "hyperbola_offset": mid + s,
        "hyperbola_offset_minus": mid - s,
    }


def symbolic():
    import sympy as sp

    c1, c2, chi, alpha = sp.symbols("c1 c2 chi alpha", real=True)
    return {
        "hyperbola_offset": (c1**2 - c2**2) / 3 + sp.sqrt(alpha**2 - chi**2) / 6,
        "hyperbola_offset_minus": (c1**2 - c2**2) / 3 - sp.sqrt(alpha**2 - chi**2) / 6,
    }
