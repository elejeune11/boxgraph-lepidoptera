import numpy as np


def _gl_integrate(func, a, b, n_quad):
    if b <= a:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    xs = mid + half * nodes
    return float(half * np.sum(weights * func(xs)))


def run(inputs, settings):
    v = float(inputs["v"])
    x_c = float(inputs["x_c"])
    L_D = float(inputs["L_D"])
    L_y = float(inputs["L_y"])
    L_M = float(inputs["L_M"])
    n_quad = max(2, int(round(float(settings.get("n_quad", 16)))))

    y_c = v * L_y
    a = x_c - 0.5 * L_M
    b = x_c + 0.5 * L_M

    def integrand(x):
        diff = L_D ** 2 - x ** 2 - y_c ** 2
        return (2.0 * np.abs(x) * abs(y_c)) / np.sqrt(4.0 * (y_c ** 2) * (L_D ** 2) + diff ** 2)

    if a < 0.0 < b:
        integral_um = _gl_integrate(integrand, a, 0.0, n_quad) + _gl_integrate(integrand, 0.0, b, n_quad)
    else:
        integral_um = _gl_integrate(integrand, a, b, n_quad)

    return {"Delta_F_over_f0": integral_um}
