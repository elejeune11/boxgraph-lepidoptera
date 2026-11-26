"""Least-squares fit of I_bg(q) = A / q to a digitized SAXS background curve."""

import numpy as np


def run(inputs, settings):
    q = np.asarray(inputs["I_curve"]["x"], dtype=float)
    y = np.asarray(inputs["I_curve"]["y"], dtype=float)
    inv_q = 1.0 / q
    A = float(np.dot(y, inv_q) / np.dot(inv_q, inv_q))
    y_fit = A * inv_q
    return {"I_fit": {"x": [float(v) for v in q], "y": [float(v) for v in y_fit]}}
