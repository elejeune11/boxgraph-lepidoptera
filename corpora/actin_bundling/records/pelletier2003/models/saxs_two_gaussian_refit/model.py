"""Nonlinear least-squares two-Gaussian fit to background-subtracted SAXS data I_sub(q)."""

import math
import numpy as np
from scipy.optimize import curve_fit


def _two_gauss(q, a1, m1, s1, a2, m2, s2):
    return a1 * np.exp(-0.5 * ((q - m1) / s1) ** 2) + a2 * np.exp(-0.5 * ((q - m2) / s2) ** 2)


def run(inputs, settings):
    q = np.asarray(inputs["I_sub"]["x"], dtype=float)
    y = np.asarray(inputs["I_sub"]["y"], dtype=float)
    p0 = [float(y.max()) * 0.8, 0.19, 0.028, float(y.max()) * 0.55, 0.29, 0.065]
    bounds = ([0.0, 0.14, 0.005, 0.0, 0.23, 0.010], [5000.0, 0.23, 0.12, 5000.0, 0.36, 0.25])
    popt, _ = curve_fit(_two_gauss, q, y, p0=p0, bounds=bounds, maxfev=10000)
    a1, m1, s1, a2, m2, s2 = [float(v) for v in popt]
    y_fit = _two_gauss(q, a1, m1, s1, a2, m2, s2)
    k_fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0))
    return {
        "I_fit": {"x": [float(v) for v in q], "y": [float(v) for v in y_fit]},
        "q_1": m1,
        "q_2": m2,
        "w_1": 2.0 * s1,
        "w_2": 2.0 * s2,
        "fwhm_1": k_fwhm * s1,
        "fwhm_2": k_fwhm * s2,
    }
