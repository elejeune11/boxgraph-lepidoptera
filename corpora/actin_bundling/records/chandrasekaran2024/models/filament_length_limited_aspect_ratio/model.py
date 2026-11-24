def run(inputs, settings):
    return {"aspect_ratio": (inputs["L_max"] / (2.0 * inputs["R_drop"])) ** 1.5}


def symbolic():
    import sympy as sp

    L_max, R_drop = sp.symbols("L_max R_drop", positive=True)
    return {"aspect_ratio": (L_max / (2 * R_drop)) ** sp.Rational(3, 2)}
