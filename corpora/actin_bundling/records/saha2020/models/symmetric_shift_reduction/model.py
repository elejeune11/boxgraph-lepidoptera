def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    return {
        "c1": abs(c12 - c11) / 2.0,
        "c2": abs(c22 - c21) / 2.0,
        "phi1_mid": (c11 + c12) / 2.0,
        "phi2_mid": (c21 + c22) / 2.0,
    }


def symbolic():
    import sympy as sp

    c11, c12, c21, c22 = sp.symbols("c11 c12 c21 c22", real=True)
    return {
        "c1": sp.Abs(c12 - c11) / 2,
        "c2": sp.Abs(c22 - c21) / 2,
        "phi1_mid": (c11 + c12) / 2,
        "phi2_mid": (c21 + c22) / 2,
    }
