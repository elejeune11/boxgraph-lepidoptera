import numpy as np


def _rhs(x1, x2, c1, c2, chi, chip, alpha):
    m1 = 4.0 * x1 * (x1 * x1 - c1 * c1) + (chi + alpha) * x2 + 2.0 * chip * x1 * x2 * x2
    m2 = 4.0 * x2 * (x2 * x2 - c2 * c2) + (chi - alpha) * x1 + 2.0 * chip * x2 * x1 * x1
    return -m1, -m2


def run(inputs, settings):
    c1 = abs(float(inputs["c12"]) - float(inputs["c11"])) / 2.0
    c2 = abs(float(inputs["c22"]) - float(inputs["c21"])) / 2.0
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])
    tg = np.array([float(v) for v in inputs["target_green"]["x"]])
    yg = np.array([float(v) for v in inputs["target_green"]["y"]])
    tb = np.array([float(v) for v in inputs["target_blue"]["x"]])
    yb = np.array([float(v) for v in inputs["target_blue"]["y"]])

    dt = float(settings.get("dt", 0.05))
    T = float(max(tg.max(), tb.max())) + 120.0
    x1 = float(settings.get("x1_init", -1e-4))
    x2 = float(settings.get("x2_init", -1e-4))

    n = int(round(T / dt)) + 1
    ts = np.arange(n) * dt
    xs1 = np.empty(n)
    xs2 = np.empty(n)
    xs1[0], xs2[0] = x1, x2
    for k in range(1, n):
        k1 = _rhs(x1, x2, c1, c2, chi, chip, alpha)
        k2 = _rhs(x1 + 0.5 * dt * k1[0], x2 + 0.5 * dt * k1[1], c1, c2, chi, chip, alpha)
        k3 = _rhs(x1 + 0.5 * dt * k2[0], x2 + 0.5 * dt * k2[1], c1, c2, chi, chip, alpha)
        k4 = _rhs(x1 + dt * k3[0], x2 + dt * k3[1], c1, c2, chi, chip, alpha)
        x1 += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        x2 += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        xs1[k], xs2[k] = x1, x2

    s_sum = xs1 + xs2
    s_diff = xs1 - xs2

    def _fit_t0(grid):
        best = None
        for t0 in grid:
            p1 = np.interp(tg + t0, ts, s_sum)
            p2 = np.interp(tb + t0, ts, s_diff)
            err = float(np.mean((p1 - yg) ** 2 + (p2 - yb) ** 2))
            if best is None or err < best[0]:
                best = (err, float(t0), p1, p2)
        return best

    b0 = _fit_t0(np.linspace(-30.0, 90.0, 121))
    b1 = _fit_t0(np.linspace(b0[1] - 1.5, b0[1] + 1.5, 41))
    return {
        "green_curve": {"x": tg.tolist(), "y": b1[2].tolist()},
        "blue_curve": {"x": tb.tolist(), "y": b1[3].tolist()},
    }
