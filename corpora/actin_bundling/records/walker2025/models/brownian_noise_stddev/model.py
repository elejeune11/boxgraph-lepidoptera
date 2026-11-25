import math


def run(inputs, settings):
    # D in um^2/s, dt in s -> sigma in um -> convert to nm (*1000)
    return {"sigma_B": math.sqrt(2.0 * float(inputs["D"]) * float(inputs["dt"])) * 1000.0}


def symbolic():
    import sympy as sp

    D, dt = sp.symbols("D dt", positive=True)
    return {"sigma_B": 1000 * sp.sqrt(2 * D * dt)}
