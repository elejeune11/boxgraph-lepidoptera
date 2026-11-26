import math


def run(inputs, settings):
    L = float(inputs["L_mem_0"])
    c_0 = float(inputs["c_0"])
    dt_total = float(inputs["dt_total"])
    f_mem = float(inputs["f_mem"])      # pN
    delta = float(inputs["delta"])      # nm
    N = float(inputs["N"])              # dimensionless
    kBT = float(inputs["kBT"])          # pN*nm
    n_steps = int(settings.get("n_steps", 200))

    rate = c_0 * math.exp(-f_mem * delta / (N * kBT))
    h = dt_total / n_steps
    for _ in range(n_steps):
        k1 = rate / (1.0 + L)
        k2 = rate / (1.0 + L + 0.5 * h * k1)
        k3 = rate / (1.0 + L + 0.5 * h * k2)
        k4 = rate / (1.0 + L + h * k3)
        L += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return {"L_mem": L}
