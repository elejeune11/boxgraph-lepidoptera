"""Single-Gaussian nonlinear least-squares fit to a digitized plotted Gaussian component curve."""

import math
import numpy as np
from scipy.optimize import curve_fit


def _one_gauss(q, A, mu, sigma):
    return A * np.exp(-0.5 * ((q - mu) / sigma) ** 2)


def run(inputs, settings):
    q = np.asarray(inputs["G_curve"]["x"], dtype=float)
    y = np.asarray(inputs["G_curve"]["y"], dtype=float)
    i_max = int(np.argmax(y))
    p0 = [float(y[i_max]), float(q[i_max]), max(0.015, 0.25 * float(q.max() - q.min()))]
    bounds = ([0.0, 0.10, 0.005], [5000.0, 0.40, 0.25])
    popt, _ = curve_fit(_one_gauss, q, y, p0=p0, bounds=bounds, maxfev=10000)
    A, mu, sigma = [float(v) for v in popt]
    y_fit = _one_gauss(q, A, mu, sigma)
    k_fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0))
    return {
        "G_fit": {"x": [float(v) for v in q], "y": [float(v) for v in y_fit]},
        "q_0": mu,
        "w_2sigma": 2.0 * sigma,
        "fwhm": k_fwhm * sigma,
    }
