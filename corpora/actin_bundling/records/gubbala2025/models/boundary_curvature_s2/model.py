import numpy as np


def run(inputs, settings):
    R = float(inputs["R"])
    n_pts = int(settings.get("n_pts", 360))
    u = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False)
    # Circle g(u) = (R cos u, R sin u)
    xp = -R * np.sin(u)
    yp = R * np.cos(u)
    xpp = -R * np.cos(u)
    ypp = -R * np.sin(u)
    num = np.abs(xp * ypp - yp * xpp)
    # Printed Eq. (S2) denominator: (x' + y')^(3/2) is real and finite only where x' + y' > 0
    # (on the other 50% of the circle it is complex/singular). Evaluate printed formula over its real domain.
    base = xp + yp
    real_mask = base > 1e-12 * R
    kappa_real = num[real_mask] / (base[real_mask] ** 1.5)
    kappa_max = float(np.max(kappa_real))
    return {
        "kappa_mean": float(np.mean(kappa_real)),
        "kappa_norm_mean": float(np.mean(kappa_real / kappa_max)),
    }
