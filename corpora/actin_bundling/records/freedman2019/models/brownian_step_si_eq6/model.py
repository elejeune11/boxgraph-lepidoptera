import math
import numpy as np

KB_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    F = float(inputs["F"])
    mu = float(inputs["mu"])
    T = float(inputs["T"])
    t_total = float(inputs["t_total"])

    dt = float(settings.get("dt", 0.00002))
    n_traj = int(settings.get("n_traj", 2000))
    seed = int(settings.get("seed", 0))

    n_steps = max(1, int(round(t_total / dt)))
    actual_t = n_steps * dt
    kBT = KB_PN_UM_PER_K * T
    sigma = math.sqrt(0.5 * kBT * mu * dt)
    drift = F * mu * dt

    rng = np.random.default_rng(seed)
    x = np.zeros(n_traj, dtype=np.float64)
    w_prev = rng.standard_normal(n_traj)
    for _ in range(n_steps):
        w_next = rng.standard_normal(n_traj)
        x += drift + sigma * (w_next + w_prev)
        w_prev = w_next

    mean_drift = F * mu * actual_t
    msd = float(np.mean((x - mean_drift) ** 2))
    return {"D_eff": msd / (2.0 * actual_t)}
