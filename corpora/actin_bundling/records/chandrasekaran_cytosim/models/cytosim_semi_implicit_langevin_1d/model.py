"""Faithful Python port of CytoSim semi-implicit Meca::solve() vs Couple::diffuse():
- src/sim/meca.cc:345-346 (tau_ = time_step; alpha_ = 2 * kT / tau_)
- src/sim/meca.cc:552-556 ((I - dt*mob*M)*dx = dt*mob*(F + Noise))
- src/sim/bead.cc:109-119 (b = sqrt(alpha * drag), rhs = fce + b * Z)
- src/sim/couple_prop.cc:146, src/sim/couple.cc:138 (dx = sqrt(6*D*dt) * sreal())
"""
import math
import numpy as np


def run(inputs, settings):
    mu = float(inputs["mu"])
    kT = float(inputs["kT"])
    k_trap = float(inputs["k_trap"])
    d_couple = float(inputs["D_couple"])

    dt = float(settings.get("dt", 0.002))
    t_sim = float(settings.get("t_sim", 1.0))
    n_particles = int(settings.get("n_particles", 4000))
    seed = int(settings.get("seed", 42))

    n_steps = max(1, int(round(t_sim / dt)))
    t_actual = n_steps * dt
    rng = np.random.default_rng(seed)

    # 1. Free Bead in Meca::solve() (M = 0, F = 0)
    # Bead::addBrownianForces: b = sqrt((2*kT/dt) / mu); rhs = dt * mu * b * Z
    alpha = 2.0 * kT / dt
    b_force = math.sqrt(alpha / mu)
    step_sigma = dt * mu * b_force  # == sqrt(2 * kT * mu * dt)
    x_free = np.zeros(n_particles, dtype=np.float64)

    # 2. Harmonically trapped Bead in Meca::solve() (M = -k_trap, F = -k_trap * x)
    # Semi-implicit factor: 1 / (1 + dt * mu * k_trap)
    inv_lhs = 1.0 / (1.0 + dt * mu * k_trap)
    x_trap = np.zeros(n_particles, dtype=np.float64)

    # 3. Free Couple in Couple::diffuse() (uniform [-sqrt(6*D*dt), +sqrt(6*D*dt)])
    hw_couple = math.sqrt(6.0 * d_couple * dt)
    x_couple = np.zeros(n_particles, dtype=np.float64)

    burn_in = n_steps // 2
    var_trap_accum = 0.0
    var_trap_count = 0

    for step in range(n_steps):
        z1 = rng.standard_normal(n_particles)
        x_free += step_sigma * z1

        z2 = rng.standard_normal(n_particles)
        rhs_trap = dt * mu * (-k_trap * x_trap + b_force * z2)
        x_trap += inv_lhs * rhs_trap
        if step >= burn_in:
            var_trap_accum += float(np.mean(x_trap ** 2))
            var_trap_count += 1

        u = rng.uniform(-1.0, 1.0, size=n_particles)
        x_couple += hw_couple * u

    d_eff_free = float(np.mean(x_free ** 2)) / (2.0 * t_actual)
    d_eff_couple = float(np.mean(x_couple ** 2)) / (2.0 * t_actual)
    var_trap = var_trap_accum / max(1, var_trap_count)
    var_boltz = kT / k_trap

    return {
        "D_eff_meca_free": d_eff_free,
        "D_eff_couple_uniform": d_eff_couple,
        "var_x_semi_implicit": var_trap,
        "var_x_boltzmann": var_boltz,
    }
