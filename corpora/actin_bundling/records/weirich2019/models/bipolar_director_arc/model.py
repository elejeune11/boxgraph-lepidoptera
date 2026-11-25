import math


def run(inputs, settings):
    x = inputs["x"]
    y = inputs["y"]
    L_D = inputs["L_D"]
    diff = L_D ** 2 - x ** 2 - y ** 2
    R = math.sqrt(L_D ** 2 + (diff ** 2) / (4.0 * y ** 2))
    sin_alpha_d = (2.0 * x * y) / math.sqrt(4.0 * (y ** 2) * (L_D ** 2) + diff ** 2)
    return {"R": R, "sin_alpha_d": sin_alpha_d}


def symbolic():
    import sympy as sp
    x, y, L_D = sp.symbols("x y L_D", positive=True)
    diff = L_D ** 2 - x ** 2 - y ** 2
    R = sp.sqrt(L_D ** 2 + diff ** 2 / (4 * y ** 2))
    sin_alpha_d = (2 * x * y) / sp.sqrt(4 * y ** 2 * L_D ** 2 + diff ** 2)
    return {"R": R, "sin_alpha_d": sin_alpha_d}
