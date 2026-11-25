import math
import numpy as np


def run(inputs, settings):
    mu = inputs["mu"]
    kBT_pN_um = inputs["kBT"] * 1e-3
    f = inputs["f"]
    dt = float(settings.get("dt", 0.001))
    t_sim = float(settings.get("t_sim", 0.2))
    n_particles = int(settings.get("n_particles", 2000))
    seed = int(settings.get("seed", 42))

    D = mu * kBT_pN_um
    n_steps = max(1, int(round(t_sim / dt)))
    actual_t = n_steps * dt
    drift_step = mu * f * dt
    sigma = math.sqrt(2.0 * D * dt)

    rng = np.random.default_rng(seed)
    x = np.zeros(n_particles, dtype=np.float64)
    for _ in range(n_steps):
        x += drift_step + rng.normal(0.0, sigma, size=n_particles)

    mean_dx = float(np.mean(x))
    var_dx = float(np.var(x, ddof=1))
    D_eff = var_dx / (2.0 * actual_t)
    return {"D_eff": D_eff, "mean_dx": mean_dx}
