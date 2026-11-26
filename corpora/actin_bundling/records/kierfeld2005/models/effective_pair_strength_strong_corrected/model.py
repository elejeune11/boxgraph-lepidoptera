import math


def run(inputs, settings):
    X_1 = inputs["X_1"]
    W = inputs["W"]
    a_par = inputs["a_par"]
    T = inputs["T"]
    z = (1.0 - X_1) + X_1 * math.exp(-a_par * W / T)
    return {"W_bar": -(T / a_par) * math.log(z)}


def symbolic():
    import sympy as sp
    X_1, W, a_par, T = sp.symbols("X_1 W a_par T", real=True)
    return {"W_bar": -(T / a_par) * sp.log((1 - X_1) + X_1 * sp.exp(-a_par * W / T))}
