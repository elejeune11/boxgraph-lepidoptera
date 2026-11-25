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
        k1 = _rhs(x1, x2, c11, c12, c21, c22, chi, chip, alpha)
        k2 = _rhs(x1 + 0.5 * dt * k1[0], x2 + 0.5 * dt * k1[1], c11, c12, c21, c22, chi, chip, alpha)
        k3 = _rhs(x1 + 0.5 * dt * k2[0], x2 + 0.5 * dt * k2[1], c11, c12, c21, c22, chi, chip, alpha)
        k4 = _rhs(x1 + dt * k3[0], x2 + dt * k3[1], c11, c12, c21, c22, chi, chip, alpha)
        x1 += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        x2 += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        xs1[k] = x1
        xs2[k] = x2

    if mode == "traj":
        # Interpolate trajectory x2(x1) if trajectory visits target x1 range
        order = np.argsort(xs1)
        x1_s, x2_s = xs1[order], xs2[order]
        out_y = []
        for xt in xs_target:
            if x1_s[0] <= xt <= x1_s[-1]:
                out_y.append(float(np.interp(xt, x1_s, x2_s)))
            else:
                out_y.append(float("nan"))
        return {"x2_curve": {"x": xs_target, "y": out_y}}

    # Steady-state limit cycle (second half)
    h = n // 2
    s1, s2 = xs1[h:], xs2[h:]
    if (s1.max() - s1.min()) < 1e-3:
        return {"x2_curve": {"x": xs_target, "y": [float("nan")] * len(xs_target)}}
    out_y = []
    for xt in xs_target:
        # Find crossings of s1 == xt
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
