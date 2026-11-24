import math
import numpy as np


def run(inputs, settings):
    mu = inputs["mu"]  # um / (pN * s)
    kbt_pnum = inputs["kBT"] * 1e-3  # pN * nm -> pN * um
    f = inputs["f"]  # pN
    dt = float(settings["dt"])
    t_sim = float(settings["t_sim"])
    n_particles = int(round(settings["n_particles"]))
    seed = int(settings["seed"])

    n_steps = max(1, int(round(t_sim / dt)))
    t_actual = n_steps * dt
    d_true = mu * kbt_pnum
    sigma = math.sqrt(2.0 * d_true * dt)
    drift = mu * f * dt

    rng = np.random.default_rng(seed)
    x = np.zeros(n_particles, dtype=np.float64)
    for _ in range(n_steps):
        x += drift + sigma * rng.standard_normal(n_particles)

    mean_dx = float(np.mean(x))
    var_dx = float(np.var(x, ddof=1))
    d_eff = var_dx / (2.0 * t_actual)
    return {"D_eff": d_eff, "mean_dx": mean_dx}
