"""Deterministic fourth-order Runge-Kutta solver of Eqs. S1, S3, S5, S8, S9 (printed sign -delta_{i1}).

Simulates two parallel actin filaments (each 2 particles) bridged by one alpha-actinin
crosslinker (1 center particle + 2 ABDs at coordinate alpha) relaxing from an initial
transverse crosslinker perturbation delta_y0 under zero thermal noise.
"""

import math
import numpy as np

CAP_UM = 100.0


def _rhs(state, params, plus_sign):
    # state: shape (5, 2) -> [f0_0, f0_1, f1_0, f1_1, cl_c]
    r00, r01, r10, r11, rc = state[0], state[1], state[2], state[3], state[4]
    L0_f = params["L0_f"]
    kappa_f = params["kappa_f"]
    L0_cl = params["L0_cl"]
    kappa_cl = params["kappa_cl"]
    k_bend = params["k_cl_bend"]
    alpha = params["alpha"]
    gx_f = params["gx_f"]
    gy_f = params["gy_f"]
    g_cl = params["g_cl"]

    forces = np.zeros_like(state)

    # Filament stretching forces (Eq. S5)
    for idx0, idx1 in ((0, 1), (2, 3)):
        d = state[idx1] - state[idx0]
        dist = math.hypot(d[0], d[1])
        if dist < 1e-15:
            return None
        f1 = -(kappa_f / L0_f) * (dist - L0_f) * (d / dist)
        forces[idx1] += f1
        forces[idx0] -= f1

    # Crosslinker ABD positions
    r_abd0 = (1.0 - alpha) * r00 + alpha * r01
    r_abd1 = (1.0 - alpha) * r10 + alpha * r11

    # Crosslinker arm vectors l_{c,h} = r_{c,h} - r_c (Eqs. S8, S9)
    l0 = r_abd0 - rc
    l1 = r_abd1 - rc
    d0 = math.hypot(l0[0], l0[1])
    d1 = math.hypot(l1[0], l1[1])
    if d0 < 1e-15 or d1 < 1e-15:
        return None

    # Arm stretching forces on ABDs (Eq. S8, rest length L0_cl/2)
    half_L = 0.5 * L0_cl
    f_abd0_str = -(kappa_cl / half_L) * (d0 - half_L) * (l0 / d0)
    f_abd1_str = -(kappa_cl / half_L) * (d1 - half_L) * (l1 / d1)

    # Bending forces on ABDs (Eq. S9): u = -l0, v = l1
    u = -l0
    v = l1
    uhat = u / d0
    vhat = v / d1
    cos_t = uhat[0] * vhat[0] + uhat[1] * vhat[1]
    f_abd0_bend = -k_bend * (vhat - cos_t * uhat) / d0
    f_abd1_bend = k_bend * (uhat - cos_t * vhat) / d1

    f_abd0 = f_abd0_str + f_abd0_bend
    f_abd1 = f_abd1_str + f_abd1_bend

    # Reaction force on crosslinker center particle (Eq. S3)
    forces[4] = -(f_abd0 + f_abd1)

    # Distribute ABD forces to filament particles (SI Notes p. 2)
    s1 = 1.0 if plus_sign else -1.0
    forces[0] += (1.0 - alpha) * f_abd0
    forces[1] += s1 * alpha * f_abd0
    forces[2] += (1.0 - alpha) * f_abd1
    forces[3] += s1 * alpha * f_abd1

    vel = np.empty_like(state)
    vel[0:4, 0] = forces[0:4, 0] / gx_f
    vel[0:4, 1] = forces[0:4, 1] / gy_f
    vel[4, :] = forces[4, :] / g_cl
    return vel


def solve_dimer(inputs, settings, plus_sign):
    dt = float(settings["dt"])
    eta = inputs["eta_cyto"]
    d_f = inputs["d_f"]
    L0_f = inputs["L0_f"]
    d_cl = inputs["d_cl"]
    L0_cl = inputs["L0_cl"]
    k_bend = inputs["k_cl_bend"]
    asym = inputs["friction_asymmetry"]
    alpha = inputs["alpha"]
    dy0 = inputs["delta_y0"]
    t_max = inputs["t_max"]

    gx_f = 3.0 * math.pi * eta * (3.0 * d_f + 0.4 * L0_f)
    gy_f = asym * gx_f
    g_cl = 6.0 * math.pi * eta * (3.0 * d_cl + 0.2 * L0_cl)

    # Dimensionless linear bending stiffness number lambda_perp * dt
    lambda_perp = 16.0 * k_bend / (L0_cl ** 2 * g_cl)
    stiffness_num = lambda_perp * dt

    params = {
        "L0_f": L0_f,
        "kappa_f": inputs["kappa_f"],
        "L0_cl": L0_cl,
        "kappa_cl": inputs["kappa_cl"],
        "k_cl_bend": k_bend,
        "alpha": alpha,
        "gx_f": gx_f,
        "gy_f": gy_f,
        "g_cl": g_cl,
    }

    a = 0.5 * L0_cl
    y_abd_init = (alpha - 0.5) * L0_f
    state = np.array(
        [
            [-a, -0.5 * L0_f],
            [-a, +0.5 * L0_f],
            [+a, -0.5 * L0_f],
            [+a, +0.5 * L0_f],
            [0.0, y_abd_init + dy0],
        ],
        dtype=float,
    )

    def weighted_cm(st):
        wx = gx_f * np.sum(st[0:4, 0]) + g_cl * st[4, 0]
        wy = gy_f * np.sum(st[0:4, 1]) + g_cl * st[4, 1]
        return np.array([wx / (4.0 * gx_f + g_cl), wy / (4.0 * gy_f + g_cl)])

    cm0 = weighted_cm(state)
    n_steps = max(1, int(round(t_max / dt)))
    h = t_max / n_steps
    diverged = False

    for _ in range(n_steps):
        k1 = _rhs(state, params, plus_sign)
        if k1 is None or np.max(np.abs(state + 0.5 * h * k1)) > CAP_UM:
            diverged = True
            break
        k2 = _rhs(state + 0.5 * h * k1, params, plus_sign)
        if k2 is None or np.max(np.abs(state + 0.5 * h * k2)) > CAP_UM:
            diverged = True
            break
        k3 = _rhs(state + 0.5 * h * k2, params, plus_sign)
        if k3 is None or np.max(np.abs(state + h * k3)) > CAP_UM:
            diverged = True
            break
        k4 = _rhs(state + h * k3, params, plus_sign)
        if k4 is None:
            diverged = True
            break
        state = state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if not np.all(np.isfinite(state)) or np.max(np.abs(state)) > CAP_UM:
            diverged = True
            break

    if diverged:
        return {
            "filament_y_shift": CAP_UM * 1000.0,
            "cm_drift": CAP_UM * 1000.0,
            "rk4_stiffness_number": stiffness_num,
        }

    r_abd0 = (1.0 - alpha) * state[0] + alpha * state[1]
    r_abd1 = (1.0 - alpha) * state[2] + alpha * state[3]
    y_shift_nm = (0.5 * (r_abd0[1] + r_abd1[1]) - y_abd_init) * 1000.0
    cm_end = weighted_cm(state)
    cm_drift_nm = math.hypot(cm_end[0] - cm0[0], cm_end[1] - cm0[1]) * 1000.0
    return {
        "filament_y_shift": float(y_shift_nm),
        "cm_drift": float(cm_drift_nm),
        "rk4_stiffness_number": float(stiffness_num),
    }


def run(inputs, settings):
    return solve_dimer(inputs, settings, plus_sign=False)
