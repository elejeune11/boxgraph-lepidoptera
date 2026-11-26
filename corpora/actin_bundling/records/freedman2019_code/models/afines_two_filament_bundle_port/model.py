"""Faithful short-horizon Python port of AFINES 2-filament bundle dynamics loop:
- prog/bundles.cpp:410-485 (main integration loop)
- src/filament.cpp:206-242, 518-575 (Leimkuhler-Matthews BD + LAMMPS harmonic bending)
- src/spring.cpp:77-81, 124-129 (Hookean link stretching)
- src/spacer.cpp:213-298, 329-337 (doubly-bound spacer tension + angular bending + lever rule)
- src/motor_ensemble.cpp:313-317 (kb = s_bend / mld)
- src/filament_ensemble.cpp:571-715 (excluded volume)

Integrates 2 parallel 3-link filaments initialized at y = ±0.125 um bridged by 9 static seed spacers
evenly spaced on link 0 (matching grow.cfg/nogrow.cfg + row_of_spacers.py geometry) over a short
mechanical relaxation window t_sim_s (default 0.01 s = 500 steps at dt = 2e-5 s) so validation and
dt-refinement complete in seconds while exercising every C++ force kernel together.
"""

import math
import numpy as np


def run(inputs, settings):
    l0 = float(inputs["link_length"])
    kl = float(inputs["link_stretching_stiffness"])
    kappa_b = float(inputs["polymer_bending_modulus"])
    rad = float(inputs["actin_length"])
    mld = float(inputs["spacer_length"])
    mk = float(inputs["spacer_stiffness"])
    s_bend = float(inputs["spacer_bend"])
    kexv = float(inputs["kexv"])
    rmax = float(inputs["rmax"])
    kT = float(inputs["temperature"])
    visc = float(inputs["viscosity"])
    dt_base = float(inputs["dt"])
    nfixed = int(round(float(inputs["nfixed"])))

    seed = int(settings.get("seed", 1))
    t_sim = float(settings.get("t_sim_s", 0.01))
    n_links = int(settings.get("n_links", 3))
    dt_scale = float(settings.get("dt_scale", 1.0))

    dt = dt_base * dt_scale
    n_steps = max(1, int(round(t_sim / dt)))
    n_beads = n_links + 1

    rng = np.random.default_rng(seed)

    # Derived C++ constants
    kb_fil = kappa_b / l0  # prog/bundles.cpp:224
    kb_sp = s_bend / mld   # src/motor_ensemble.cpp:314
    damp = 6.0 * math.pi * rad * visc  # src/bead.cpp:25
    bd_pref = math.sqrt(kT / (2.0 * dt * damp))  # src/filament.cpp:83
    max_small = 0.001  # include/globals.h:58
    th0 = 0.5 * math.pi  # prog/bundles.cpp:123

    # Initialize 2 parallel filaments along -x starting at x=-11.5, y=±0.125 (grow.cfg:12)
    pos = np.zeros((2, n_beads, 2), dtype=float)
    for f in range(2):
        y_init = -0.125 if f == 0 else 0.125
        for b in range(n_beads):
            pos[f, b, 0] = -11.5 - b * l0
            pos[f, b, 1] = y_init
    prv_rnd = np.zeros_like(pos)

    # 9 seed spacers evenly spaced along link 0: pos_a_end = distance from bead 1 (pointed end of link 0)
    # row_of_spacers.py places 9 spacers at x = -11.5 - 0.1*k (k=0..8), so distance from bead 1 (-12.5) is 1.0 - 0.1*k
    # clamp slightly inside link to avoid r1=0 singularity at bead 0
    pos_a = np.array([max(0.05, min(0.95, l0 - 0.1 * (k + 0.5))) for k in range(nfixed)], dtype=float)

    pe_tot = 0.0
    for _ in range(n_steps):
        forces = np.zeros_like(pos)
        pe_stretch = 0.0
        pe_bend = 0.0
        pe_sp = 0.0
        pe_exv = 0.0

        # 1. Filament stretching & LAMMPS bending
        for f in range(2):
            disp = pos[f, 1:] - pos[f, :-1]  # shape (n_links, 2)
            lens = np.hypot(disp[:, 0], disp[:, 1])
            direc = disp / lens[:, None]
            for l in range(n_links):
                kf = kl * (lens[l] - l0)
                fx, fy = kf * direc[l, 0], kf * direc[l, 1]
                forces[f, l, 0] += fx
                forces[f, l, 1] += fy
                forces[f, l + 1, 0] -= fx
                forces[f, l + 1, 1] -= fy
                pe_stretch += 0.5 * kl * ((lens[l] - l0) ** 2)

            for n in range(n_links - 1):
                dr1 = -disp[n]
                dr2 = disp[n + 1]
                r1 = lens[n]
                r2 = lens[n + 1]
                c = max(-1.0, min(1.0, (dr1[0] * dr2[0] + dr1[1] * dr2[1]) / (r1 * r2)))
                s = math.sqrt(max(0.0, 1.0 - c * c))
                if s < max_small:
                    s = max_small
                theta = math.acos(c) - math.pi
                pe_bend += 0.5 * kb_fil * theta * theta
                a = -kb_fil * theta / s
                a11 = a * c / (r1 * r1)
                a12 = -a / (r1 * r2)
                a22 = a * c / (r2 * r2)
                f1x = a11 * dr1[0] + a12 * dr2[0]
                f1y = a11 * dr1[1] + a12 * dr2[1]
                f3x = a22 * dr2[0] + a12 * dr1[0]
                f3y = a22 * dr2[1] + a12 * dr1[1]
                forces[f, n, 0] += f1x
                forces[f, n, 1] += f1y
                forces[f, n + 1, 0] += -f1x - f3x
                forces[f, n + 1, 1] += -f1y - f3y
                forces[f, n + 2, 0] += f3x
                forces[f, n + 2, 1] += f3y

        # 2. Doubly-bound seed spacers on link 0 (src/spacer.cpp + src/motor.cpp:486-491)
        d0 = pos[0, 1] - pos[0, 0]
        L0 = math.hypot(d0[0], d0[1])
        u0 = d0 / L0
        d1 = pos[1, 1] - pos[1, 0]
        L1 = math.hypot(d1[0], d1[1])
        u1 = d1 / L1

        for pa in pos_a:
            h0 = pos[0, 1] - pa * u0
            h1 = pos[1, 1] - pa * u1
            sdisp = h1 - h0
            slen = math.hypot(sdisp[0], sdisp[1])
            if slen <= 1e-12:
                continue
            sdir = sdisp / slen
            tens = mk * (slen - mld)
            ft = tens * sdir
            pe_sp += 0.5 * mk * ((slen - mld) ** 2)

            # Bending at head 0 and head 1
            bf0 = np.zeros(2)
            bf1 = np.zeros(2)
            for hd in (0, 1):
                h_cur = h0 if hd == 0 else h1
                L_cur = L0 if hd == 0 else L1
                fe = 1 if pa > 0.5 * L_cur else 0
                dr1 = pos[hd, fe] - h_cur
                r1 = math.hypot(dr1[0], dr1[1])
                dr2 = sdisp if hd == 0 else -sdisp
                r2 = slen
                c = max(-1.0, min(1.0, (dr1[0] * dr2[0] + dr1[1] * dr2[1]) / (r1 * r2)))
                s = math.sqrt(max(0.0, 1.0 - c * c))
                if s < max_small:
                    s = max_small
                dth = math.acos(c) - th0
                pe_sp += 0.5 * kb_sp * dth * dth
                a = -kb_sp * dth / s
                a11 = a * c / (r1 * r1)
                a12 = -a / (r1 * r2)
                a22 = a * c / (r2 * r2)
                f1v = np.array([a11 * dr1[0] + a12 * dr2[0], a11 * dr1[1] + a12 * dr2[1]])
                f3v = np.array([a22 * dr2[0] + a12 * dr1[0], a22 * dr2[1] + a12 * dr1[1]])
                forces[hd, fe] += f1v
                if hd == 0:
                    bf0 += -f1v - f3v
                    bf1 += f3v
                else:
                    bf1 += -f1v - f3v
                    bf0 += f3v

            r_rat0 = pa / L0
            fh0 = ft + bf0
            forces[0, 0] += fh0 * r_rat0
            forces[0, 1] += fh0 * (1.0 - r_rat0)

            r_rat1 = pa / L1
            fh1 = -ft + bf1
            forces[1, 0] += fh1 * r_rat1
            forces[1, 1] += fh1 * (1.0 - r_rat1)

        # 3. Excluded volume between beads across filaments when kexv > 0 (closest-point check)
        if kexv > 0.0:
            inv_rmax = 1.0 / rmax
            for b in range(n_beads):
                dv = pos[1, b] - pos[0, b]
                r = math.hypot(dv[0], dv[1])
                if 1e-12 < r < rmax:
                    pe_exv += kexv * ((1.0 - r * inv_rmax) ** 2)
                    fac = 2.0 * kexv * inv_rmax * ((1.0 / r) - inv_rmax)
                    forces[1, b] += fac * dv
                    forces[0, b] -= fac * dv

        # 4. Leimkuhler-Matthews overdamped BD position update (src/filament.cpp:219-236)
        new_rnd = rng.standard_normal(size=pos.shape)
        vel = forces / damp + bd_pref * (new_rnd + prv_rnd)
        pos = pos + vel * dt
        prv_rnd = new_rnd
        pe_tot = pe_stretch + pe_bend + pe_sp + pe_exv

    e2e0 = pos[0, -1] - pos[0, 0]
    e2e1 = pos[1, -1] - pos[1, 0]
    cos_th = float(np.dot(e2e0, e2e1) / (np.linalg.norm(e2e0) * np.linalg.norm(e2e1)))
    mean_sep = float(np.mean(np.hypot(pos[1, :, 0] - pos[0, :, 0], pos[1, :, 1] - pos[0, :, 1])))

    return {
        "final_alignment": cos_th,
        "final_mean_separation": mean_sep,
        "final_total_pe": float(pe_tot),
    }
