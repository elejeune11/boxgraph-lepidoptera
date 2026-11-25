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
    xs_target = [float(v) for v in inputs["target"]["x"]]

    mode = str(settings.get("mode", "lc_upper"))
    dt = float(settings.get("dt", 0.02))
    T = float(settings.get("T", 200.0))
    x1 = float(settings.get("x1_init", -0.001))
    x2 = float(settings.get("x2_init", -0.001))

    n = int(round(T / dt))
    xs1 = np.empty(n)
    xs2 = np.empty(n)
    for k in range(n):
        k1 = _rhs(x1, x2, c1, c2, chi, chip, alpha)
        k2 = _rhs(x1 + 0.5 * dt * k1[0], x2 + 0.5 * dt * k1[1], c1, c2, chi, chip, alpha)
        k3 = _rhs(x1 + 0.5 * dt * k2[0], x2 + 0.5 * dt * k2[1], c1, c2, chi, chip, alpha)
        k4 = _rhs(x1 + dt * k3[0], x2 + dt * k3[1], c1, c2, chi, chip, alpha)
        x1 += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        x2 += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        xs1[k] = x1
        xs2[k] = x2

    if mode == "traj":
        # Monotonic trajectory from near origin into lower-left node
        order = np.argsort(xs1)
        x1_s, x2_s = xs1[order], xs2[order]
        out_y = [float(np.interp(xt, x1_s, x2_s)) if x1_s[0] <= xt <= x1_s[-1] else float("nan") for xt in xs_target]
        return {"x2_curve": {"x": xs_target, "y": out_y}}

    h = n // 2
    s1, s2 = xs1[h:], xs2[h:]
    if (s1.max() - s1.min()) < 1e-3:
        return {"x2_curve": {"x": xs_target, "y": [float("nan")] * len(xs_target)}}
    out_y = []
    for xt in xs_target:
        cands = []
        for i in range(len(s1) - 1):
            a, b = s1[i] - xt, s1[i + 1] - xt
            if (a <= 0.0 < b) or (b <= 0.0 < a):
                f = -a / (b - a)
                cands.append(float(s2[i] + f * (s2[i + 1] - s2[i])))
        if not cands:
            out_y.append(float("nan"))
        else:
            out_y.append(max(cands) if mode == "lc_upper" else min(cands))
    return {"x2_curve": {"x": xs_target, "y": out_y}}
