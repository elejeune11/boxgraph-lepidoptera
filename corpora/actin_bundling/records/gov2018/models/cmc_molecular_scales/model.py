def run(inputs, settings):
    abs_H0 = abs(inputs["H0"])
    abs_H = abs(inputs["H"])
    a_um = 1.0 / abs_H0
    R_um = 1.0 / abs_H
    return {
        "a": a_um * 1000.0,
        "ns": abs_H0 ** 2,
        "R": R_um,
        "a_over_R": abs_H / abs_H0,
    }


def symbolic():
    import sympy as sp

    H0, H = sp.symbols("H0 H", positive=True)
    return {
        "a": 1000 / H0,
        "ns": H0 ** 2,
        "R": 1 / H,
        "a_over_R": H / H0,
    }
