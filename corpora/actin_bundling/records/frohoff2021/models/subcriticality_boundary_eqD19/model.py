import math
from scipy.optimize import brentq


def run(inputs, settings):
    Sigma = inputs["Sigma"]
    M = abs(inputs["M"])
    Delta_H = 0.25 * M * M

    if Sigma <= 0.0 or M <= 0.0:
        return {"Delta_sub_plus": None, "Delta_sub_minus": None, "Delta_H": Delta_H}

    # Blue region (n+ subcriticality): exists for Sigma > Delta_H; lower boundary solves xi_+(Delta)^2 / Delta = Sigma
    Delta_plus = None
    if Sigma > Delta_H:
        def f_plus(d):
            xi = (0.5 * M + math.sqrt(max(0.0, Delta_H - d))) ** 2
            return xi * xi / d - Sigma
        Delta_plus = float(brentq(f_plus, 1e-12, Delta_H, xtol=1e-12, rtol=1e-12))

    # Orange region (n- subcriticality): exists for 0 < Sigma < Delta_H; lower/right boundary solves xi_-(Delta)^2 / Delta = Sigma
    Delta_minus = None
    if 0.0 < Sigma < Delta_H:
        def f_minus(d):
            xi = (0.5 * M - math.sqrt(max(0.0, Delta_H - d))) ** 2
            return xi * xi / d - Sigma
        Delta_minus = float(brentq(f_minus, 1e-12, Delta_H, xtol=1e-12, rtol=1e-12))

    return {
        "Delta_sub_plus": Delta_plus,
        "Delta_sub_minus": Delta_minus,
        "Delta_H": Delta_H,
    }
