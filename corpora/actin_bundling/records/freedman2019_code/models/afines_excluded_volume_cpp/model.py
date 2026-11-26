"""Faithful line-by-line Python port of AFINES excluded volume in src/filament_ensemble.cpp:571-715
and segment geometry helpers in src/spring.cpp:261-358.

Demonstrates:
1. For non-intersecting parallel links separated by r < rmax, Fx1 = 2*kexv*dist*(1/rmax)*(1/r - 1/rmax)
   is the exact negative gradient of pe_exv = kexv*(1 - r/rmax)^2, proving kexv has units pN*um (energy),
   NOT pN*um^2 as misprinted in freedman2019 Table S1 and Fig. 2D.
2. When two 1-um links cross orthogonally at their midpoints (r_c[0..3] = L/2 > rmax), the outer
   `if (r < rmax)` check skips the block completely, giving ZERO excluded-volume force.
3. When two links intersect with r < rmax, the `else` branch (lines 700-713) applies a hardcoded
   +45-degree force Fx1 = Fy1 = 2*kexv/(rmax*sqrt(2)) to BOTH beads of link 1 (total magnitude 4*kexv/rmax).
"""

import math


def _get_r_c(a0, a1, px, py):
    """Port of spring::get_r_c (src/spring.cpp:261-312) in non-wrapping local coordinates."""
    dx_s = a1[0] - a0[0]
    dy_s = a1[1] - a0[1]
    l2 = dx_s * dx_s + dy_s * dy_s
    if l2 == 0.0:
        pt = a0
    else:
        tp = ((px - a0[0]) * dx_s + (py - a0[1]) * dy_s) / l2
        if tp < 0.0:
            pt = a0
        elif tp > 1.0:
            pt = a1
        else:
            pt = (a0[0] + tp * dx_s, a0[1] + tp * dy_s)
    return math.hypot(px - pt[0], py - pt[1]), pt


def _get_line_intersect(a0, a1, b0, b1):
    """Port of spring::get_line_intersect (src/spring.cpp:314-358)."""
    disp1 = (a1[0] - a0[0], a1[1] - a0[1])
    disp2 = (b1[0] - b0[0], b1[1] - b0[1])
    disp12 = (a0[0] - b0[0], a0[1] - b0[1])
    denom = disp1[0] * disp2[1] - disp1[1] * disp2[0]
    if denom == 0.0:
        return False
    denom_pos = denom > 0.0
    s_num = disp1[0] * disp12[1] - disp1[1] * disp12[0]
    t_num = disp2[0] * disp12[1] - disp2[1] * disp12[0]
    if (s_num < 0.0) == denom_pos:
        return False
    if (t_num < 0.0) == denom_pos:
        return False
    if ((s_num > denom) == denom_pos) or ((t_num > denom) == denom_pos):
        return False
    return True


def _eval_exv_pair(a0, a1, b0, b1, kexv, rmax):
    """Port of filament_ensemble::update_force_between_filaments (src/filament_ensemble.cpp:571-715)."""
    b = 1.0 / rmax
    r0, p1 = _get_r_c(a0, a1, b0[0], b0[1])
    r1, p2 = _get_r_c(a0, a1, b1[0], b1[1])
    r2, p3 = _get_r_c(b0, b1, a0[0], a0[1])
    r3, p4 = _get_r_c(b0, b1, a1[0], a1[1])
    rc = [r0, r1, r2, r3]
    r = rc[0]
    idx = 0
    for k in range(1, 4):
        if rc[k] < r:
            r = rc[k]
            idx = k
    intersect = _get_line_intersect(a0, a1, b0, b1)

    len_a = math.hypot(a1[0] - a0[0], a1[1] - a0[1])
    len_b = math.hypot(b1[0] - b0[0], b1[1] - b0[1])
    f_a0 = [0.0, 0.0]
    f_a1 = [0.0, 0.0]
    f_b0 = [0.0, 0.0]
    f_b1 = [0.0, 0.0]
    pe = 0.0

    if r < rmax:
        if not intersect:
            if idx == 0:
                x1, y1, x2, y2 = b0[0], b0[1], p1[0], p1[1]
                len1, length = math.hypot(a0[0] - x2, a0[1] - y2), len_a
            elif idx == 1:
                x1, y1, x2, y2 = b1[0], b1[1], p2[0], p2[1]
                len1, length = math.hypot(a0[0] - x2, a0[1] - y2), len_a
            elif idx == 2:
                x1, y1, x2, y2 = a0[0], a0[1], p3[0], p3[1]
                len1, length = math.hypot(b0[0] - x2, b0[1] - y2), len_b
            else:
                x1, y1, x2, y2 = a1[0], a1[1], p4[0], p4[1]
                len1, length = math.hypot(b0[0] - x2, b0[1] - y2), len_b

            dist_x, dist_y = x2 - x1, y2 - y1
            len2 = length - len1
            w1 = len2 / length
            w2 = len1 / length
            fx1 = 2.0 * kexv * dist_x * b * ((1.0 / r) - b)
            fy1 = 2.0 * kexv * dist_y * b * ((1.0 / r) - b)
            pe = kexv * ((1.0 - r * b) ** 2)
            if idx in (0, 1):
                f_a0 = [fx1 * w1, fy1 * w1]
                f_a1 = [fx1 * w2, fy1 * w2]
                if idx == 0:
                    f_b0 = [-fx1, -fy1]
                else:
                    f_b1 = [-fx1, -fy1]
            else:
                f_b0 = [fx1 * w1, fy1 * w1]
                f_b1 = [fx1 * w2, fy1 * w2]
                if idx == 2:
                    f_a0 = [-fx1, -fy1]
                else:
                    f_a1 = [-fx1, -fy1]
        else:
            fx1 = 2.0 * kexv / (rmax * math.sqrt(2.0))
            fy1 = 2.0 * kexv / (rmax * math.sqrt(2.0))
            pe = kexv * ((1.0 - r * b) ** 2)
            f_a0 = [fx1, fy1]
            f_a1 = [fx1, fy1]
            f_b0 = [-fx1, -fy1]
            f_b1 = [-fx1, -fy1]

    f_tot_a = math.hypot(f_a0[0] + f_a1[0], f_a0[1] + f_a1[1])
    return pe, f_tot_a


def run(inputs, settings):
    kexv = float(inputs["kexv"])
    rmax = float(inputs["rmax"])
    sep = float(inputs["sep_parallel"])
    L = float(inputs["link_length"])

    # Case 1: parallel links separated by sep (< rmax)
    u_par, f_par = _eval_exv_pair((0.0, 0.0), (L, 0.0), (0.0, sep), (L, sep), kexv, rmax)
    f_analytic = 2.0 * kexv / rmax * (1.0 - sep / rmax) if sep < rmax else 0.0

    # Case 2: orthogonal crossing at midpoints (r = L/2 > rmax -> zero force bug!)
    _, f_cross_mid = _eval_exv_pair((-0.5 * L, 0.0), (0.5 * L, 0.0), (0.0, -0.5 * L), (0.0, 0.5 * L), kexv, rmax)

    # Case 3: orthogonal crossing near endpoint (distance sep < rmax from endpoint -> hardcoded 45-deg doubled force!)
    _, f_cross_end = _eval_exv_pair((0.0, 0.0), (L, 0.0), (sep, -0.5 * L), (sep, 0.5 * L), kexv, rmax)

    return {
        "U_parallel": u_par,
        "F_parallel_code": f_par,
        "F_parallel_analytic_neg_grad": f_analytic,
        "F_crossing_midpoint_code": f_cross_mid,
        "F_crossing_near_end_total_code": f_cross_end,
    }
