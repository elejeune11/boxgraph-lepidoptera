"""Stochastic 2D overdamped Langevin solver for a single 2-particle actin filament (Eqs. S1a-S1b, S5).

Uses Euler-Maruyama with anisotropic friction (gamma_x^f = gamma^f, gamma_y^f = asym * gamma^f)
and fluctuation-dissipation theorem <xi_{i,a}(0) xi_{j,b}(t)> = 2 gamma_a^f k_B T delta_{ij} delta_{ab} delta(t).
"""

import math
import numpy as np


def run(inputs, settings):
    dt = float(settings["dt"])
    n_traj = int(settings["n_traj"])
    seed = int(settings["seed"])
    rng = np.random.default_rng(seed)

    eta = inputs["eta_cyto"]
    d_f = inputs["d_f"]
    L0_f = inputs["L0_f"]
    kappa_f = inputs["kappa_f"]
    kBT = inputs["kBT"]
    asym = inputs["friction_asymmetry"]
    t_obs = inputs["t_obs"]

    gx = 3.0 * math.pi * eta * (3.0 * d_f + 0.4 * L0_f)
    gy = asym * gx
    k_eff = kappa_f / L0_f

    n_steps = max(1, int(round(t_obs / dt)))
    h = t_obs / n_steps

    # Thermal kick standard deviations on particle positions: sqrt(2 kBT h / gamma_a)
    sx = math.sqrt(2.0 * kBT * h / gx)
    sy = math.sqrt(2.0 * kBT * h / gy)

    # Random initial filament orientations centered at origin
    theta0 = rng.uniform(0.0, 2.0 * math.pi, size=n_traj)
    half = 0.5 * L0_f
    r0 = np.column_stack((-half * np.cos(theta0), -half * np.sin(theta0)))
    r1 = np.column_stack((+half * np.cos(theta0), +half * np.sin(theta0)))

    for _ in range(n_steps):
        d = r1 - r0
        dist = np.hypot(d[:, 0], d[:, 1])
        dist_safe = np.maximum(dist, 1e-12)
        f_mag = -k_eff * (dist_safe - L0_f) / dist_safe
        fx = f_mag * d[:, 0]
        fy = f_mag * d[:, 1]

        nx0 = rng.standard_normal(n_traj) * sx
        ny0 = rng.standard_normal(n_traj) * sy
        nx1 = rng.standard_normal(n_traj) * sx
        ny1 = rng.standard_normal(n_traj) * sy

        r1[:, 0] += (fx / gx) * h + nx1
        r1[:, 1] += (fy / gy) * h + ny1
        r0[:, 0] -= (fx / gx) * h - nx0
        r0[:, 1] -= (fy / gy) * h - ny0

    cm = 0.5 * (r0 + r1)
    mean_dx_nm = float(np.mean(np.abs(cm[:, 0])) * 1000.0)
    mean_dy_nm = float(np.mean(np.abs(cm[:, 1])) * 1000.0)
    return {
        "mean_abs_dx": mean_dx_nm,
        "mean_abs_dy": mean_dy_nm,
        "R_aniso": float(mean_dy_nm / mean_dx_nm),
    }
