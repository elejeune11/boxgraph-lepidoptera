"""Faithful line-by-line Python port of AFINES spacer bending, force propagation, and Metropolis MC:
- src/motor_ensemble.cpp:313-317 (kb = s_bend / mld)
- src/spacer.cpp:213-298, 329-337 (spacer::update_force, update_bending, filament_update)
- src/spacer.cpp:346-382 (spacer::metropolis_prob)
- src/motor.cpp:486-491 (motor::filament_update_hd lever rule)

Quantifies three implementation bugs in src/spacer.cpp:220-224, 355-378:
1. When head 1 attaches while head 0 is bound, spacer::metropolis_prob() evaluates delr1 and delr2
   from pre-attachment (hx[1], hy[1]) and old disp rather than newpos, and queries stale pos_a_end[1]=0.
2. It adds bend_eng for only ONE head (head 1), omitting the new bending energy simultaneously
   activated at head 0 upon entering state (1, 1).
3. spacer::update_force() lines 220-224 never zero b_eng when unbinding, so stale b_eng[hd] from a
   prior doubly-bound state is erroneously subtracted in delE line 377.
"""

import math


def _bend_angle_energy(delr1, delr2, kb, th0):
    r1 = math.hypot(delr1[0], delr1[1])
    r2 = math.hypot(delr2[0], delr2[1])
    c = (delr1[0] * delr2[0] + delr1[1] * delr2[1]) / (r1 * r2)
    c = max(-1.0, min(1.0, c))
    dth = math.acos(c) - th0
    return 0.5 * kb * dth * dth, c, r1, r2, dth


def _lammps_angle_forces(delr1, delr2, kb, th0, max_small_angle=0.001):
    eng, c, r1, r2, dth = _bend_angle_energy(delr1, delr2, kb, th0)
    s = math.sqrt(max(0.0, 1.0 - c * c))
    if s < max_small_angle:
        s = max_small_angle
    a = -kb * dth / s
    a11 = a * c / (r1 * r1)
    a12 = -a / (r1 * r2)
    a22 = a * c / (r2 * r2)
    f1 = (a11 * delr1[0] + a12 * delr2[0], a11 * delr1[1] + a12 * delr2[1])
    f3 = (a22 * delr2[0] + a12 * delr1[0], a22 * delr2[1] + a12 * delr1[1])
    return eng, f1, f3


def run(inputs, settings):
    mld = float(inputs["mld"])
    mk = float(inputs["mk"])
    s_bend = float(inputs["s_bend"])
    th0 = float(inputs["th0"])
    L = float(inputs["link_length"])
    sep_y = float(inputs["head_sep_y"])
    tilt_dx = float(inputs["tilt_dx"])
    pre_dy = float(inputs["pre_bind_dy"])
    pos_a_end = float(inputs["pos_a_end"])
    stale_b_eng = float(inputs["stale_b_eng"])
    kon = float(inputs["kon"])
    dt = float(inputs["dt"])
    kT = float(inputs["temperature"])

    # Effective angular bending stiffness (src/motor_ensemble.cpp:314)
    kb = s_bend / mld

    # Two parallel horizontal links along +x:
    # Filament 0: bead 0 at (0, 0), bead 1 at (L, 0). Head 0 bound at x0 = L - pos_a_end, y0 = 0.
    # Filament 1: bead 0 at (0, sep_y), bead 1 at (L, sep_y).
    # Candidate attachment point on filament 1: newpos = (x0 + tilt_dx, sep_y)
    # Pre-attachment unbound head 1 position: old_h1 = (x0 + tilt_dx, sep_y - pre_dy)
    x0 = L - pos_a_end
    h0 = (x0, 0.0)
    newpos = (x0 + tilt_dx, sep_y)
    old_h1 = (x0 + tilt_dx, sep_y - pre_dy)

    f0_beads = [(0.0, 0.0), (L, 0.0)]
    f1_beads = [(0.0, sep_y), (L, sep_y)]

    # --- A. Doubly-bound equilibrium forces & conservation check (spacer::update_force + filament_update) ---
    disp_bound = (newpos[0] - h0[0], newpos[1] - h0[1])
    len_bound = math.hypot(disp_bound[0], disp_bound[1])
    tension = mk * (len_bound - mld)
    f_tens = (tension * disp_bound[0] / len_bound, tension * disp_bound[1] / len_bound)

    pos_a_0 = L - h0[0]
    pos_a_1 = L - newpos[0]
    fe0 = 1 if pos_a_0 > 0.5 * L else 0
    fe1 = 1 if pos_a_1 > 0.5 * L else 0
    delr1_0 = (f0_beads[fe0][0] - h0[0], f0_beads[fe0][1] - h0[1])
    delr2_0 = (disp_bound[0], disp_bound[1])
    eng0, f1_0, f3_0 = _lammps_angle_forces(delr1_0, delr2_0, kb, th0)

    delr1_1 = (f1_beads[fe1][0] - newpos[0], f1_beads[fe1][1] - newpos[1])
    delr2_1 = (-disp_bound[0], -disp_bound[1])
    eng1, f1_1, f3_1 = _lammps_angle_forces(delr1_1, delr2_1, kb, th0)

    b_force0 = (-f1_0[0] - f3_0[0] + f3_1[0], -f1_0[1] - f3_0[1] + f3_1[1])
    b_force1 = (f3_0[0] - f1_1[0] - f3_1[0], f3_0[1] - f1_1[1] - f3_1[1])

    # Accumulate all forces applied to the 4 beads (f0_0, f0_1, f1_0, f1_1)
    bead_f = [[0.0, 0.0] for _ in range(4)]
    bead_f[fe0][0] += f1_0[0]
    bead_f[fe0][1] += f1_0[1]
    bead_f[2 + fe1][0] += f1_1[0]
    bead_f[2 + fe1][1] += f1_1[1]

    r0 = pos_a_0 / L
    fh0 = (f_tens[0] + b_force0[0], f_tens[1] + b_force0[1])
    bead_f[0][0] += fh0[0] * r0
    bead_f[0][1] += fh0[1] * r0
    bead_f[1][0] += fh0[0] * (1.0 - r0)
    bead_f[1][1] += fh0[1] * (1.0 - r0)

    r1 = pos_a_1 / L
    fh1 = (-f_tens[0] + b_force1[0], -f_tens[1] + b_force1[1])
    bead_f[2][0] += fh1[0] * r1
    bead_f[2][1] += fh1[1] * r1
    bead_f[3][0] += fh1[0] * (1.0 - r1)
    bead_f[3][1] += fh1[1] * (1.0 - r1)

    fx_sum = sum(bf[0] for bf in bead_f)
    fy_sum = sum(bf[1] for bf in bead_f)
    f_res = math.hypot(fx_sum, fy_sum)

    # --- B. Second-head attachment Metropolis probability: C++ buggy vs true physical Delta U ---
    disp_old = (old_h1[0] - h0[0], old_h1[1] - h0[1])
    len_old = math.hypot(disp_old[0], disp_old[1])
    u_stretch_old = 0.5 * mk * ((len_old - mld) ** 2)
    u_stretch_new = 0.5 * mk * ((len_bound - mld) ** 2)

    # As written in src/spacer.cpp:355-378:
    # - get_further_end(1, ...) checks pos_a_end[1] (= 0 while unbound!) > 0.5*L -> False -> bead 0!
    # - disp_from_bead(1, 1, 0) uses pre-attachment old_h1 instead of newpos!
    # - delr2 uses pre-attachment -disp_old instead of -disp_bound!
    # - subtracts unreset stale_b_eng for head 1, and omits head 0's newly activated eng0!
    delr1_bug = (f1_beads[0][0] - old_h1[0], f1_beads[0][1] - old_h1[1])
    delr2_bug = (-disp_old[0], -disp_old[1])
    bend_eng_bug, _, _, _, _ = _bend_angle_energy(delr1_bug, delr2_bug, kb, th0)
    delE_bug = u_stretch_new + bend_eng_bug - u_stretch_old - stale_b_eng

    # True physical mechanical energy change (state (1,0) has 0 bending energy; state (1,1) has eng0 + eng1)
    delE_true = (u_stretch_new - u_stretch_old) + (eng0 + eng1)

    maxprob = kon * dt
    prob_bug = maxprob * math.exp(-delE_bug / kT) if delE_bug > 0.0 else maxprob
    prob_true = maxprob * math.exp(-delE_true / kT) if delE_true > 0.0 else maxprob

    return {
        "kb_eff": kb,
        "U_bend_both_heads": eng0 + eng1,
        "F_total_residual": f_res,
        "delE_code_buggy": delE_bug,
        "delE_true_physical": delE_true,
        "delE_discrepancy": delE_true - delE_bug,
        "prob_attach_code_buggy": prob_bug,
        "prob_attach_true_physical": prob_true,
    }
