"""Faithful line-by-line Python port of AFINES filament mechanics:
- src/spring.cpp:77-81, 376-378 (Hookean stretching force & energy)
- src/filament.cpp:518-575 (LAMMPS harmonic angle bending force & energy)
- src/filament.cpp:80-83, 206-242 (Leimkuhler-Matthews overdamped BD step)
"""

import math


def run(inputs, settings):
    kl = float(inputs["kl"])
    l0 = float(inputs["l0"])
    llen = float(inputs["llen"])
    kappa_b = float(inputs["polymer_bending_modulus"])
    ang = float(inputs["bend_angle"])
    max_small_angle = float(inputs["max_small_angle"])
    kT = float(inputs["temperature"])
    visc = float(inputs["viscosity"])
    rad = float(inputs["actin_length"])
    dt = float(inputs["dt"])

    # 1. Hookean link stretching (src/spring.cpp:79-80, 377)
    f_stretch = kl * (llen - l0)
    u_stretch = (f_stretch * f_stretch) / (2.0 * kl)

    # 2. 3-bead geometry: bead 0 = (-l0, 0), bead 1 = (0, 0), bead 2 = (l0*cos(ang), l0*sin(ang))
    #    spring 0: bead 0 -> bead 1, disp = (l0, 0); delr1 = neg_disp = (-l0, 0)
    #    spring 1: bead 1 -> bead 2, disp = (l0*cos(ang), l0*sin(ang))
    kb = kappa_b / l0  # prog/bundles.cpp:224
    delr1 = (-l0, 0.0)
    delr2 = (l0 * math.cos(ang), l0 * math.sin(ang))
    rsq1 = delr1[0] ** 2 + delr1[1] ** 2
    rsq2 = delr2[0] ** 2 + delr2[1] ** 2
    r1 = math.sqrt(rsq1)
    r2 = math.sqrt(rsq2)
    c = (delr1[0] * delr2[0] + delr1[1] * delr2[1]) / (r1 * r2)
    c = max(-1.0, min(1.0, c))
    s = math.sqrt(max(0.0, 1.0 - c * c))
    if s < max_small_angle:
        s = max_small_angle
    theta = math.acos(c) - math.pi
    u_bend = 0.5 * kb * theta * theta

    a = -kb * theta / s
    a11 = a * c / rsq1
    a12 = -a / (r1 * r2)
    a22 = a * c / rsq2
    f1 = (a11 * delr1[0] + a12 * delr2[0], a11 * delr1[1] + a12 * delr2[1])
    f3 = (a22 * delr2[0] + a12 * delr1[0], a22 * delr2[1] + a12 * delr1[1])
    f2 = (-f1[0] - f3[0], -f1[1] - f3[1])

    f_bend_end = math.hypot(f1[0], f1[1])
    f_net_res = math.hypot(f1[0] + f2[0] + f3[0], f1[1] + f2[1] + f3[1])

    # 3. Leimkuhler-Matthews effective 1D diffusion coefficient (src/filament.cpp:80, 83, 222-233)
    #    dx_brn = bd_prefactor * dt * (W_n + W_{n-1}), bd_prefactor = sqrt(kT / (2*dt*damp))
    #    Long-time variance per coordinate per dt is 2 * (bd_prefactor * dt)^2 * 2 = 2 * (kT / damp) * dt = 2 * D * dt
    damp = 6.0 * math.pi * rad * visc
    bd_prefactor = math.sqrt(kT / (2.0 * dt * damp))
    d_eff = 2.0 * (bd_prefactor * bd_prefactor) * dt  # equals kT / damp exactly

    return {
        "U_stretch": u_stretch,
        "F_stretch": f_stretch,
        "U_bend": u_bend,
        "F_bend_end_bead": f_bend_end,
        "F_bend_net_residual": f_net_res,
        "D_bead_lm": d_eff,
    }
