"""Line-by-line Python port of AFINES Leimkuhler-Matthews (2013) overdamped Langevin update:
- src/bead.cpp:20-22 / src/motor.cpp:35 (damp = 6 * pi * vis * r)
- src/filament.cpp:67 / src/motor.cpp:36 (bd_prefactor = sqrt(temperature / (2 * damp * dt)))
- src/motor.cpp:97-98 (prv_rnd initialized to N(0, 1) to prevent cooling)
- src/filament.cpp:95-101 / src/motor.cpp:485-491 (v = F/damp + bd_prefactor*(new_rnd + prv_rnd); pos += v*dt)
"""

import math
import numpy as np


def run(inputs, settings):
    vis = float(inputs["viscosity"])  # mg/(um*s) = pN*s/um^2
    r_hydro = float(inputs["r_hydro"])
    temperature = float(inputs["temperature"])
    k_trap = float(inputs["k_trap"])
    dt = float(inputs["dt"])

    n_steps = int(round(float(settings.get("n_steps", 2000))))
    n_particles = int(round(float(settings.get("n_particles", 2000))))
    seed = int(round(float(settings.get("seed", 42))))

    damp = 6.0 * math.pi * vis * r_hydro
    bd_prefactor = math.sqrt(temperature / (2.0 * damp * dt))
    D_theory = temperature / damp
    var_x_boltzmann = temperature / k_trap

    rng = np.random.default_rng(seed)

    # 1. Free 2D diffusion with exact AFINES Leimkuhler-Matthews update
    pos_free = np.zeros((n_particles, 2), dtype=float)
    prv_rnd_free = rng.standard_normal((n_particles, 2))
    for _ in range(n_steps):
        new_rnd = rng.standard_normal((n_particles, 2))
        v = bd_prefactor * (new_rnd + prv_rnd_free)
        pos_free += v * dt
        prv_rnd_free = new_rnd
    msd_2d = float(np.mean(pos_free[:, 0] ** 2 + pos_free[:, 1] ** 2))
    t_total = n_steps * dt
    D_sim = msd_2d / (4.0 * t_total)

    # 2. Harmonically trapped 2D particle (F = -k_trap * pos, matching ext_circle excess coordinate)
    pos_trap = rng.normal(0.0, math.sqrt(var_x_boltzmann), size=(n_particles, 2))
    prv_rnd_trap = rng.standard_normal((n_particles, 2))
    burn = n_steps // 2
    var_acc = 0.0
    var_cnt = 0
    for step in range(n_steps):
        force = -k_trap * pos_trap
        new_rnd = rng.standard_normal((n_particles, 2))
        v = force / damp + bd_prefactor * (new_rnd + prv_rnd_trap)
        pos_trap += v * dt
        prv_rnd_trap = new_rnd
        if step >= burn:
            var_acc += float(np.mean(pos_trap[:, 0] ** 2))
            var_cnt += 1
    var_x_sim = var_acc / max(1, var_cnt)

    return {
        "D_theory": D_theory,
        "D_sim": D_sim,
        "var_x_boltzmann": var_x_boltzmann,
        "var_x_sim": var_x_sim,
    }
