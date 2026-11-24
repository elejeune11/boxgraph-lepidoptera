import math


def run(inputs, settings):
    S = float(inputs["S"])
    theta = float(inputs["theta"])
    nx = math.cos(theta)
    ny = math.sin(theta)
    Q_xx = S * (nx * nx - 0.5)
    Q_xy = S * nx * ny
    Q_yy = S * (ny * ny - 0.5)
    tr_Q2 = Q_xx * Q_xx + 2.0 * Q_xy * Q_xy + Q_yy * Q_yy
    return {"Q_xx": Q_xx, "Q_xy": Q_xy, "Q_yy": Q_yy, "tr_Q2": tr_Q2}


def symbolic():
    import sympy as sp

    S, theta = sp.symbols("S theta", real=True)
    nx = sp.cos(theta)
    ny = sp.sin(theta)
    Q_xx = sp.simplify(S * (nx**2 - sp.Rational(1, 2)))
    Q_xy = sp.simplify(S * nx * ny)
    Q_yy = sp.simplify(S * (ny**2 - sp.Rational(1, 2)))
    tr_Q2 = sp.simplify(Q_xx**2 + 2 * Q_xy**2 + Q_yy**2)
    return {"Q_xx": Q_xx, "Q_xy": Q_xy, "Q_yy": Q_yy, "tr_Q2": tr_Q2}
