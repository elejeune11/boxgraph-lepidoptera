def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    return {
        "f1_pp": a + 3.0 * p1**2,
        "f2_pp": a + a_Delta + 3.0 * p2**2,
    }


def symbolic():
    import sympy as sp

    a, a_Delta, phi1_bar, phi2_bar = sp.symbols("a a_Delta phi1_bar phi2_bar", real=True)
    return {
        "f1_pp": a + 3 * phi1_bar**2,
        "f2_pp": a + a_Delta + 3 * phi2_bar**2,
    }
