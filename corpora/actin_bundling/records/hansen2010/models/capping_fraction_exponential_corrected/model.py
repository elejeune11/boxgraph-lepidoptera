import math


def run(inputs, settings):
    f_free = math.exp(-inputs["K_obs"] * inputs["t"])
    return {"f_capped": 1.0 - f_free, "f_free": f_free}


def symbolic():
    import sympy as sp

    K_obs, t = sp.symbols("K_obs t", positive=True)
    f_free = sp.exp(-K_obs * t)
    return {"f_capped": 1 - f_free, "f_free": f_free}
