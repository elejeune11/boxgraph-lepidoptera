import math


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

    dt = float(settings.get("dt", 0.05))
    T = float(settings.get("T", 160.0))
    x1 = float(settings.get("x1_init", 0.001))
    x2 = float(settings.get("x2_init", 0.001))

    nsteps = int(round(T / dt))
    half = nsteps // 2
    xs1 = []
    xs2 = []
    ts = []

    for step in range(nsteps):
        k1_1, k1_2 = _rhs(x1, x2, c1, c2, chi, chip, alpha)
        k2_1, k2_2 = _rhs(x1 + 0.5 * dt * k1_1, x2 + 0.5 * dt * k1_2, c1, c2, chi, chip, alpha)
        k3_1, k3_2 = _rhs(x1 + 0.5 * dt * k2_1, x2 + 0.5 * dt * k2_2, c1, c2, chi, chip, alpha)
        k4_1, k4_2 = _rhs(x1 + dt * k3_1, x2 + dt * k3_2, c1, c2, chi, chip, alpha)
        x1 += (dt / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1)
        x2 += (dt / 6.0) * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2)
        if step >= half:
            xs1.append(x1)
            xs2.append(x2)
            ts.append((step + 1) * dt)

    x1_amp = 0.5 * (max(xs1) - min(xs1))
    x2_amp = 0.5 * (max(xs2) - min(xs2))
    if x1_amp < 1e-4:
        return {"omega": 0.0, "x1_amp": x1_amp, "x2_amp": x2_amp}

    mean1 = 0.5 * (max(xs1) + min(xs1))
    crossings = []
    for i in range(len(xs1) - 1):
        a = xs1[i] - mean1
        b = xs1[i + 1] - mean1
        if a <= 0.0 < b:
            frac = -a / (b - a)
            crossings.append(ts[i] + frac * (ts[i + 1] - ts[i]))
    if len(crossings) < 2:
        omega = 0.0
    else:
        period = (crossings[-1] - crossings[0]) / (len(crossings) - 1)
        omega = 2.0 * math.pi / period
    return {"omega": omega, "x1_amp": x1_amp, "x2_amp": x2_amp}
