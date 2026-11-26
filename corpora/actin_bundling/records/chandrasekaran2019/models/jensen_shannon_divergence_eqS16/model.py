import math
import numpy as np


def run(inputs, settings):
    m1, s1 = inputs["mu_1"], inputs["sigma_1"]
    m2, s2 = inputs["mu_2"], inputs["sigma_2"]
    n_grid = int(settings.get("n_grid", 1000))
    n_sigma = float(settings.get("n_sigma", 10.0))
    x_lo = min(m1 - n_sigma * s1, m2 - n_sigma * s2)
    x_hi = max(m1 + n_sigma * s1, m2 + n_sigma * s2)
    x = np.linspace(x_lo, x_hi, n_grid)
    p = np.exp(-0.5 * ((x - m1) / s1) ** 2) / (s1 * math.sqrt(2.0 * math.pi))
    q = np.exp(-0.5 * ((x - m2) / s2) ** 2) / (s2 * math.sqrt(2.0 * math.pi))
    m = 0.5 * (p + q)
    # Integrand using natural log (log base unstated in Eq. S17); guard underflow where p or q == 0.
    term_p = np.where(p > 1e-300, p * np.log(np.where(p > 1e-300, p / np.maximum(m, 1e-300), 1.0)), 0.0)
    term_q = np.where(q > 1e-300, q * np.log(np.where(q > 1e-300, q / np.maximum(m, 1e-300), 1.0)), 0.0)
    d_pm = float(np.trapezoid(term_p, x))
    d_qm = float(np.trapezoid(term_q, x))
    return {"JSD": 0.5 * d_pm + 0.5 * d_qm}
