def run(inputs, settings):
    a = 0.5 * inputs["L_max"]
    r = inputs["R_drop"]
    b = (r * r) / a
    return {"b": b, "aspect_ratio": a / b}


def symbolic():
    import sympy as sp

    L_max, R_drop = sp.symbols("L_max R_drop", positive=True)
    a = L_max / 2
    b = R_drop**2 / a
    return {"b": b, "aspect_ratio": (L_max / (2 * R_drop)) ** 2}
