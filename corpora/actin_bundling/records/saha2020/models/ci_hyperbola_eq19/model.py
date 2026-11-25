def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])
    alpha = float(inputs["alpha"])
    return {"hyperbola_offset": (c1 * c1 + c2 * c2 - chi * chi + alpha * alpha) / 3.0}


def symbolic():
    import sympy as sp

    c1, c2, chi, alpha = sp.symbols("c1 c2 chi alpha", real=True)
    return {"hyperbola_offset": (c1**2 + c2**2 - chi**2 + alpha**2) / 3}
