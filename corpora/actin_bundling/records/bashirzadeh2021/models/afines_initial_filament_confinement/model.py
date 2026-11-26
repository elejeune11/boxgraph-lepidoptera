import math
import numpy as np


def run(inputs, settings):
    n_beads = int(round(float(inputs["N_beads"])))
    l_link = float(inputs["l_link"])
    sigma = float(inputs["sigma_disp"])
    r_c = float(inputs["r_c"])
    k_c = float(inputs["k_c"])

    n_samples = int(round(settings.get("n_samples", 500)))
    seed = int(round(settings.get("seed", 42)))
    rng = np.random.default_rng(seed)

    half_n = (n_beads - 1) // 2
    other_n = (n_beads - 1) - half_n
    arm1_k = np.arange(1, half_n + 1, dtype=float) * l_link
    arm2_k = np.arange(1, other_n + 1, dtype=float) * l_link

    theta1 = rng.uniform(0.0, 2.0 * math.pi, size=n_samples)
    delta = rng.uniform(-math.pi / 4.0, math.pi / 4.0, size=n_samples)
    theta2 = theta1 + math.pi + delta
    dx = rng.normal(0.0, sigma, size=n_samples)
    dy = rng.normal(0.0, sigma, size=n_samples)

    # Center bead + arm1 + arm2 positions for all n_samples filaments
    x0 = dx[:, None]
    y0 = dy[:, None]
    x1 = dx[:, None] + np.cos(theta1)[:, None] * arm1_k[None, :]
    y1 = dy[:, None] + np.sin(theta1)[:, None] * arm1_k[None, :]
    x2 = dx[:, None] + np.cos(theta2)[:, None] * arm2_k[None, :]
    y2 = dy[:, None] + np.sin(theta2)[:, None] * arm2_k[None, :]

    xs = np.concatenate([x0, x1, x2], axis=1)
    ys = np.concatenate([y0, y1, y2], axis=1)
    rs = np.sqrt(xs ** 2 + ys ** 2)

    excess = np.maximum(rs - r_c, 0.0)
    u_per_filament = 0.5 * k_c * np.sum(excess ** 2, axis=1)
    frac_outside = float(np.mean(rs >= r_c))
    mean_u = float(np.mean(u_per_filament))
    return {"frac_outside": frac_outside, "mean_U_confine": mean_u}
