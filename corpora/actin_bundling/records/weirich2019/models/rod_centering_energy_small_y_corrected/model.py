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
    f_0 = inputs["f_0"]      # pN
    x_c = inputs["x_c"]      # um
    y_c = inputs["y_c"]      # um
    L_D = inputs["L_D"]      # um
    L_M = inputs["L_M"]      # um
    n_quad = max(2, int(round(float(settings.get("n_quad", 16)))))

    a = x_c - 0.5 * L_M
    b = x_c + 0.5 * L_M

    def integrand(x):
        return 2.0 * np.abs(x) / (L_D ** 2 - x ** 2)

    if a < 0.0 < b:
        integral = _gl_integrate(integrand, a, 0.0, n_quad) + _gl_integrate(integrand, 0.0, b, n_quad)
    else:
        integral = _gl_integrate(integrand, a, b, n_quad)

    F_0_pN_nm = f_0 * abs(y_c) * 1000.0
    return {
        "Delta_F": F_0_pN_nm * integral,
        "Delta_F_over_F0": integral,
    }
