import numpy as np


def _rhs(x1, x2, c11, c12, c21, c22, chi, chip, alpha):
    m1 = 2.0 * (x1 - c11) * (x1 - c12) * (2.0 * x1 - c11 - c12) + (chi + alpha) * x2 + 2.0 * chip * x1 * x2 * x2
    m2 = 2.0 * (x2 - c21) * (x2 - c22) * (2.0 * x2 - c21 - c22) + (chi - alpha) * x1 + 2.0 * chip * x2 * x1 * x1
    return -m1, -m2


def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])
    tg = [float(v) for v in inputs["target_green"]["x"]]
    tb = [float(v) for v in inputs["target_blue"]["x"]]

    dt = float(settings.get("dt", 0.02))
    T = max(max(tg), max(tb)) + 2.0
    x1 = float(settings.get("x1_init", -0.001))
    x2 = float(settings.get("x2_init", -0.001))

    n = int(round(T / dt)) + 1
    ts = np.arange(n) * dt
    xs1 = np.empty(n)
    xs2 = np.empty(n)
    xs1[0], xs2[0] = x1, x2
    for k in range(1, n):
        k1 = _rhs(x1, x2, c11, c12, c21, c22, chi, chip, alpha)
        k2 = _rhs(x1 + 0.5 * dt * k1[0], x2 + 0.5 * dt * k1[1], c11, c12, c21, c22, chi, chip, alpha)
        k3 = _rhs(x1 + 0.5 * dt * k2[0], x2 + 0.5 * dt * k2[1], c11, c12, c21, c22, chi, chip, alpha)
        k4 = _rhs(x1 + dt * k3[0], x2 + dt * k3[1], c11, c12, c21, c22, chi, chip, alpha)
        x1 += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        x2 += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        xs1[k], xs2[k] = x1, x2

    return {
        "green_curve": {"x": tg, "y": [float(v) for v in np.interp(tg, ts, xs1)]},
        "blue_curve": {"x": tb, "y": [float(v) for v in np.interp(tb, ts, xs2)]},
    }
