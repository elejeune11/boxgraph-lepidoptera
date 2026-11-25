import numpy as np
from scipy.optimize import curve_fit


def run(inputs, settings):
    t = np.asarray(inputs["L_vs_t"]["x"], dtype=float)
    L = np.asarray(inputs["L_vs_t"]["y"], dtype=float)
    t0 = float(t[0])
    L_f0 = float(L[-1])
    A0 = max(0.1, float(L[0] - L_f0))
    tau0 = max(1e-3, 0.2 * float(t[-1] - t0))
    popt, _ = curve_fit(
        lambda x, lf, amp, tau: lf + amp * np.exp(-(x - t0) / tau),
        t,
        L,
        p0=[L_f0, A0, tau0],
        bounds=([0.0, 0.0, 1e-5], [np.inf, np.inf, np.inf]),
    )
    lf, amp, tau = float(popt[0]), float(popt[1]), float(popt[2])
    return {"tau": tau, "L_f": lf, "Delta_L": amp}
