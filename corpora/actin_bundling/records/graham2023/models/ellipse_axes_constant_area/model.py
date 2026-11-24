import math


def run(inputs, settings):
    R = 0.5 * inputs["D_droplet"]
    s = math.sqrt(inputs["AR"])
    return {"a": R * s, "b": R / s}


def symbolic():
    import sympy as sp

    D_droplet, AR = sp.symbols("D_droplet AR", positive=True)
    return {
        "a": (D_droplet / 2) * sp.sqrt(AR),
        "b": (D_droplet / 2) / sp.sqrt(AR),
    }
