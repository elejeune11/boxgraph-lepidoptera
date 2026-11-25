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
    f_0 = inputs["f_0"]      # pN (= pN*um / um)
    x_c = inputs["x_c"]      # um
    y_c = inputs["y_c"]      # um
    L_D = inputs["L_D"]      # um
    L_M = inputs["L_M"]      # um
    n_quad = max(2, int(round(float(settings.get("n_quad", 16)))))

    a = x_c - 0.5 * L_M
    b = x_c + 0.5 * L_M

    def integrand(x):
        diff = L_D ** 2 - x ** 2 - y_c ** 2
        return (2.0 * np.abs(x) * abs(y_c)) / np.sqrt(4.0 * (y_c ** 2) * (L_D ** 2) + diff ** 2)

    if a < 0.0 < b:
        integral_um = _gl_integrate(integrand, a, 0.0, n_quad) + _gl_integrate(integrand, 0.0, b, n_quad)
    else:
        integral_um = _gl_integrate(integrand, a, b, n_quad)

    # f_0 [pN] * integral [um] = pN*um = 1000 pN*nm
    delta_F_pN_nm = f_0 * integral_um * 1000.0
    delta_F_over_F0 = integral_um / abs(y_c) if abs(y_c) > 0 else 0.0
    return {
        "Delta_F": delta_F_pN_nm,
        "Delta_F_over_F0": delta_F_over_F0,
    }
