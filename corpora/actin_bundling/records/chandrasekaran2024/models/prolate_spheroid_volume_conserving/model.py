import math


def run(inputs, settings):
    a = inputs["a"]
    r = inputs["R_drop"]
    b = math.sqrt((r**3) / a)
    ar = a / b
    ratio_sq = (b / a) ** 2
    if ratio_sq >= 1.0 - 1e-15:
        area = 4.0 * math.pi * r * r
    else:
        e = math.sqrt(1.0 - ratio_sq)
        area = 2.0 * math.pi * b * b * (1.0 + (a / (b * e)) * math.asin(e))
    return {"b": b, "aspect_ratio": ar, "surface_area": area}


def symbolic():
    import sympy as sp

    a, R_drop = sp.symbols("a R_drop", positive=True)
    b = sp.sqrt(R_drop**3 / a)
    ar = (a / R_drop) ** sp.Rational(3, 2)
    e = sp.sqrt(1 - (R_drop / a) ** 3)
    area = 2 * sp.pi * b**2 * (1 + (a / (b * e)) * sp.asin(e))
    return {"b": b, "aspect_ratio": ar, "surface_area": area}
