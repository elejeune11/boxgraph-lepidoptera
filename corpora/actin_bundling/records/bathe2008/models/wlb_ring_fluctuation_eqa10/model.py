"""Exact Fourier mode-sum periodic ring bending stiffness from Appendix Eq. (A10)."""

import math


def run(inputs, settings):
    kappa_f = inputs["kappa_f"]
    N = inputs["N"]
    alpha = inputs["alpha"]
    chi2 = inputs["chi2"]
    n_modes = int(settings.get("n_modes", 200))

    sqrt_N = math.sqrt(N)
    # 2 * sum_{j=1}^oo 1 / (2*pi*j)^4 = 1/720
    s = 0.0
    for j in range(1, n_modes + 1):
        qL = 2.0 * math.pi * j
        c_q = (qL * qL) / 12.0
        kb_star_q = N * (1.0 + chi2 * (N - 1.0) / (1.0 + c_q * (N + sqrt_N) / alpha))
        s += 2.0 / ((qL ** 4) * kb_star_q)

    kappa_B_star = 1.0 / (720.0 * s)
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}
