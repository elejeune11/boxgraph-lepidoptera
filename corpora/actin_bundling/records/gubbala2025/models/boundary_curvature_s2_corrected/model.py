import numpy as np


def run(inputs, settings):
    R = float(inputs["R"])
    n_pts = int(settings.get("n_pts", 360))
    u = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False)
    xp = -R * np.sin(u)
    yp = R * np.cos(u)
    xpp = -R * np.cos(u)
    ypp = -R * np.sin(u)
    num = np.abs(xp * ypp - yp * xpp)
    den = (xp * xp + yp * yp) ** 1.5
    kappa = num / den
    kappa_max = float(np.max(kappa))
    return {
        "kappa_mean": float(np.mean(kappa)),
        "kappa_norm_mean": float(np.mean(kappa / kappa_max)),
    }


def symbolic():
    import sympy as sp

    R = sp.symbols("R", positive=True)
    return {"kappa_mean": 1 / R, "kappa_norm_mean": sp.Integer(1)}
