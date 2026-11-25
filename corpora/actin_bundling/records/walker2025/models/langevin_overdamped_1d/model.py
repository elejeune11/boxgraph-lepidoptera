import math
import numpy as np


def run(inputs, settings):
    mu = float(inputs["mu"])              # um / (pN * s)
    kBT_pN_um = float(inputs["kBT"]) * 1e-3  # pN * nm -> pN * um
    f = float(inputs["f"])                # pN

    dt = float(settings.get("dt", 0.002))
    t_sim = float(settings.get("t_sim", 1.0))
    n_particles = int(settings.get("n_particles", 4000))
    seed = int(settings.get("seed", 42))

    n_steps = max(1, int(round(t_sim / dt)))
    t_actual = n_steps * dt

    D = mu * kBT_pN_um                    # um^2 / s
    sigma = math.sqrt(2.0 * D * dt)       # um
    drift = mu * f * dt                   # um

    rng = np.random.default_rng(seed)
    x = np.zeros(n_particles, dtype=np.float64)
    for _ in range(n_steps):
        x += drift + rng.normal(0.0, sigma, size=n_particles)

    mean_dx = float(np.mean(x))
    var_x = float(np.var(x, ddof=1))
    D_eff = var_x / (2.0 * t_actual)
    return {"D_eff": D_eff, "mean_dx": mean_dx}
