import numpy as np
from scipy.optimize import curve_fit


def run(inputs, settings):
    t = np.asarray(inputs["AR_vs_t"]["x"], dtype=float)
    y = np.asarray(inputs["AR_vs_t"]["y"], dtype=float)
    t0 = float(t[0])
    y_eq0 = float(y[-1])
    A0 = max(0.05, float(y[0] - y_eq0))
    tau0 = max(0.5, 0.15 * float(t[-1] - t0))
    popt, _ = curve_fit(
        lambda x, yeq, amp, tau: yeq + amp * np.exp(-(x - t0) / tau),
        t,
        y,
        p0=[y_eq0, A0, tau0],
        bounds=([0.0, 0.0, 1e-4], [2.0, 2.0, np.inf]),
    )
    return {"tau": float(popt[2]), "AR_eq": float(popt[0])}
