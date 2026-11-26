import math


def run(inputs, settings):
    # Independent-segment grand-canonical trace for N = 2 bound filaments of n_seg segments (V_r = 0, V_a = W):
    # exp(-H_bar_2 / T) = [(1 - X_1) + X_1 * exp(-a_par * W / T)]^(2 * n_seg)
    X_1 = inputs["X_1"]
    W = inputs["W"]
    a_par = inputs["a_par"]
    n_seg = inputs["n_seg"]
    T = inputs["T"]
    factor = (1.0 - X_1) + X_1 * math.exp(-a_par * W / T)
    return {"H_bar_2": -2.0 * n_seg * T * math.log(factor)}


def symbolic():
    import sympy as sp
    X_1, W, a_par, n_seg, T = sp.symbols("X_1 W a_par n_seg T", real=True)
    return {"H_bar_2": -2 * n_seg * T * sp.log((1 - X_1) + X_1 * sp.exp(-a_par * W / T))}
