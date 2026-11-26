import math


def run(inputs, settings):
    # Evaluates Eq. (4) verbatim for N = 2 bound filaments with n_seg segments inside the attractive well
    # (V_r = 0, V_a = W): the printed formula places sum_l inside the exponential without prod_l outside.
    X_1 = inputs["X_1"]
    W = inputs["W"]
    a_par = inputs["a_par"]
    n_seg = inputs["n_seg"]
    T = inputs["T"]
    factor = (1.0 - X_1) + X_1 * math.exp(-n_seg * a_par * W / T)
    return {"H_bar_2": -2.0 * T * math.log(factor)}


def symbolic():
    import sympy as sp
    X_1, W, a_par, n_seg, T = sp.symbols("X_1 W a_par n_seg T", real=True)
    return {"H_bar_2": -2 * T * sp.log((1 - X_1) + X_1 * sp.exp(-n_seg * a_par * W / T))}
