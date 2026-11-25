import math


def run(inputs, settings):
    x = inputs["x"]
    x_0 = inputs["x_0"]
    delta = inputs["delta"]
    I_0 = inputs["I_0"]
    I_B = inputs["I_B"]
    val = I_0 + (I_B / (delta * math.sqrt(2.0 * math.pi))) * math.exp(-((x - x_0) ** 2) / (2.0 * delta**2))
    return {"I": val}


def symbolic():
    import sympy as sp
    x, x_0 = sp.symbols("x x_0", real=True)
    delta, I_0, I_B = sp.symbols("delta I_0 I_B", positive=True)
    return {"I": I_0 + (I_B / (delta * sp.sqrt(2 * sp.pi))) * sp.exp(-((x - x_0) ** 2) / (2 * delta**2))}
