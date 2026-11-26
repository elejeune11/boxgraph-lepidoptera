"""Exact numerical quadrature of maxian2022 Eq. S13 Gibbs-Boltzmann equilibrium distribution
P(x) = Z x^2 exp(-K_c (x - ell)^2 / (2 k_B T)) for benchmarking the code's integrators."""

import numpy as np


def run(inputs, settings):
    K_c = float(inputs["K_c"])
    ell = float(inputs["ell"])
    kBT = float(inputs["kBT"])
    n_quad = int(settings.get("n_quad", 20001))

    sig = np.sqrt(kBT / K_c)
    xs = np.linspace(0.0, ell + 12.0 * sig, n_quad)
    unnorm = (xs ** 2) * np.exp(-K_c * ((xs - ell) ** 2) / (2.0 * kBT))
    Z = float(np.trapezoid(unnorm, xs))
    pdf = unnorm / Z
    mean_x = float(np.trapezoid(xs * pdf, xs))
    var_x = float(np.trapezoid(((xs - mean_x) ** 2) * pdf, xs))
    return {
        "mean_x": mean_x,
        "std_x": float(np.sqrt(var_x)),
    }
