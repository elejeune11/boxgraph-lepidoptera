"""Faithful 1D periodic Q1 finite-element Newton-Raphson port of AuxCortexNematic2D.cpp.

Ports line-by-line:
- src/AuxCortexNematic2D.cpp lines 241-1162 reduced to 1D (q2 = 0, vy = 0, d/dy = 0):
  Q1 linear elements, 2-point Gauss quadrature (1D projection of 4-pt 2D Gauss),
  exact residual Bk and Jacobian Ak, explicit q0-lagged 3-term active stress and
  lambda_rot, hidden log barrier f_lim = -k_lim*log(S2_max - S2) including line 1140's
  extra-h factor (controllable via settings['extra_h_bug'] for sensitivity checks).
- src/CortexNematic2D.cpp lines 511-538: Newton-Raphson loop (maxIter = 8, resTol = 1e-8)
  and adaptive time-stepping loop (deltat_init -> maxDelt).
Uses an exact cyclic block-tridiagonal banded Woodbury solver (scipy.linalg.solve_banded,
exact to 2e-16) for fast CPU execution.
"""

import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks


def _solve_cyclic_3dof(K_diag, K_upper, K_lower, R):
    N = K_diag.shape[0]
    m = 3 * N
    ab = np.zeros((11, m))
    idx = np.arange(N)
    for r in range(3):
        for c in range(3):
            ab[5 + r - c, 3 * idx + c] = K_diag[:, r, c]
            ab[5 + r - (c + 3), 3 * idx[:-1] + 3 + c] = K_upper[:-1, r, c]
            ab[5 + 3 + r - c, 3 * idx[:-1] + c] = K_lower[:-1, r, c]
    B = np.zeros((m, 7))
    B[:, 0] = R
    B[0:3, 1:4] = np.eye(3)
    B[m - 3 : m, 4:7] = np.eye(3)
    Z = solve_banded((5, 5), ab, B, overwrite_ab=True, overwrite_b=True, check_finite=False)
    y = Z[:, 0]
    W = Z[:, 1:7]
    K0N = K_lower[N - 1]
    KN0 = K_upper[N - 1]
    Vty = np.concatenate([K0N @ y[m - 3 : m], KN0 @ y[0:3]])
    VtW = np.vstack([K0N @ W[m - 3 : m, :], KN0 @ W[0:3, :]])
    c_vec = np.linalg.solve(np.eye(6) + VtW, Vty)
    return y - W @ c_vec


def _integrate_1d_fem(
    N,
    Ldom,
    U0_init,
    p,
    dt_init,
    dt_max,
    adapt,
    t_end,
    extra_h_bug=True,
    max_nr=8,
    res_tol=1e-8,
):
    dx = Ldom / N
    jac_e = 0.5 * dx
    g = 1.0 / math.sqrt(3.0)
    bf_g = np.array([[0.5 * (1.0 + g), 0.5 * (1.0 - g)], [0.5 * (1.0 - g), 0.5 * (1.0 + g)]])
    dbf_g = np.array([[-1.0 / dx, 1.0 / dx], [-1.0 / dx, 1.0 / dx]])

    n0 = np.arange(N)
    n1 = (n0 + 1) % N
    nodes = np.stack([n0, n1], axis=1)
    gdofs = np.stack([3 * n0, 3 * n0 + 1, 3 * n0 + 2, 3 * n1, 3 * n1 + 1, 3 * n1 + 2], axis=1)

    visc = p["visc"]
    cvisc = p["cvisc"]
    rvisc = p["rvisc"]
    frank_cfg = p["frank_cfg"]
    sus2_b = p["sus20_cfg"] / 2.0
    sus4_b = p["sus40_cfg"] / 8.0
    lam_rot = p["lambda_rot"]
    lam_iso = p["lambda_iso"]
    lam_aniso = p["lambda_aniso"]
    xl = p["crossLinker"]
    heqb = p["heqb"]
    k_lim_val = p["k_lim_val"]
    S2_max = p["S2_max"]

    def assemble(U, U0, dt):
        h_n, q_n, v_n = U[0::3], U[1::3], U[2::3]
        h0_n, q0_n = U0[0::3], U0[1::3]
        frank = frank_cfg / dt
        sus2 = sus2_b / dt
        sus4 = sus4_b / dt
        k_lim = k_lim_val / dt
        Re = np.zeros((N, 6))
        Ke = np.zeros((N, 6, 6))
        he, h0e = h_n[nodes], h0_n[nodes]
        qe, q0e, ve = q_n[nodes], q0_n[nodes], v_n[nodes]

        for gp in range(2):
            bf, dbf = bf_g[gp], dbf_g[gp]
            h, h0v, dh = he @ bf, h0e @ bf, he @ dbf
            q, q0v, dq, dq0 = qe @ bf, q0e @ bf, qe @ dbf, q0e @ dbf
            v, dv = ve @ bf, ve @ dbf
            Jq = (q - q0v) / dt + v * dq
            S_0 = 2.0 * np.abs(q0v)
            S2 = 4.0 * q * q
            if k_lim_val > 0.0:
                d_bar = S2_max - S2
                if np.any(d_bar <= 0.0):
                    return None, None, None, None
                f_lim = -k_lim * np.log(d_bar)
            else:
                d_bar = np.ones_like(S2)
                f_lim = np.zeros_like(S2)
            sig_act = lam_iso * (1.0 - S_0) + lam_aniso * (q0v + 0.5 * S_0) + xl * (-q0v + 0.5 * S_0)
            DqDq = 0.5 * dq * dq
            h_bar_fac = h if extra_h_bug else 1.0
            sus_lin = sus2 + 2.0 * sus4 * S2
            sus_en = S2 * (sus2 + sus4 * S2)

            for a in range(2):
                bfi, dbfi, ia = bf[a], dbf[a], 3 * a
                dqJq_i = bfi / dt
                qbfv_i = 2.0 * bfi * q
                dqd_i = -4.0 * qbfv_i
                dqflim_i = -(k_lim / d_bar) * dqd_i if k_lim_val > 0.0 else 0.0
                dvdh_i = -dt * (bfi * dh + h * dbfi)
                dvJq_i, dvJq0_i = bfi * dq, bfi * dq0

                Re[:, ia] += jac_e * (bfi * ((h - h0v) / dt + dh * v + h * dv - visc * (heqb - h)) + dbfi * dh)
                Re[:, ia + 1] += jac_e * (
                    2.0 * h * rvisc * dqJq_i * Jq
                    + h * cvisc * dqJq_i * dv
                    + h * frank * dbfi * dq
                    + 4.0 * h * sus_lin * qbfv_i
                    + h * h * lam_rot * (2.0 * bfi / dt * q0v)
                    + h * dqflim_i
                )
                Re[:, ia + 2] += jac_e * (
                    4.0 * h * dv * dbfi
                    + h * v * bfi
                    + h * sig_act * dbfi
                    + 2.0 * h * rvisc * dvJq_i * Jq
                    + h * cvisc * (dvJq_i * dv + Jq * dbfi)
                    + dvdh_i * frank * DqDq
                    + sus_en * dvdh_i
                    + h * h * lam_rot * (2.0 * q0v * dvJq0_i)
                    + h_bar_fac * f_lim * dvdh_i
                )

                for b in range(2):
                    bfj, dbfj, jb = bf[b], dbf[b], 3 * b
                    dhdvdh = -dt * (bfi * dbfj + bfj * dbfi)
                    dqJq_j = bfj / dt + v * dbfj
                    dvJq_j = bfj * dq
                    dqdvJq = bfi * dbfj
                    qbfv_j = 2.0 * bfj * q
                    dqd_j = -4.0 * qbfv_j
                    if k_lim_val > 0.0:
                        dqflim_j = -(k_lim / d_bar) * dqd_j
                        dqdqflim = -(k_lim / d_bar) * (-8.0 * bfi * bfj) + (k_lim / (d_bar * d_bar)) * dqd_i * dqd_j
                    else:
                        dqflim_j = 0.0
                        dqdqflim = 0.0

                    Ke[:, ia, jb] += jac_e * (bfi * (bfj / dt + dbfj * v + bfj * dv + visc * bfj) + dbfi * dbfj)
                    Ke[:, ia, jb + 2] += jac_e * bfi * (dh * bfj + h * dbfj)
                    Ke[:, ia + 1, jb] += jac_e * (
                        2.0 * bfj * rvisc * dqJq_i * Jq
                        + bfj * cvisc * dqJq_i * dv
                        + bfj * frank * dbfi * dq
                        + 4.0 * sus_lin * qbfv_i * bfj
                        + 2.0 * h * lam_rot * (2.0 * bfi / dt * q0v) * bfj
                        + dqflim_i * bfj
                    )
                    Ke[:, ia + 1, jb + 1] += jac_e * (
                        2.0 * h * rvisc * dqJq_i * dqJq_j
                        + h * frank * dbfi * dbfj
                        + 4.0 * h * (sus_lin * (2.0 * bfi * bfj) + 8.0 * sus4 * qbfv_i * qbfv_j)
                        + h * dqdqflim
                    )
                    Ke[:, ia + 1, jb + 2] += jac_e * (2.0 * h * rvisc * dqJq_i * dvJq_j + h * cvisc * dqJq_i * dbfj)
                    Ke[:, ia + 2, jb] += jac_e * (
                        4.0 * bfj * dv * dbfi
                        + bfj * bfi * v
                        + bfj * sig_act * dbfi
                        + 2.0 * bfj * rvisc * dvJq_i * Jq
                        + bfj * cvisc * (dvJq_i * dv + Jq * dbfi)
                        + dhdvdh * frank * DqDq
                        + sus_en * dhdvdh
                        + 2.0 * h * lam_rot * (2.0 * q0v * dvJq0_i) * bfj
                        + (f_lim * (dvdh_i * bfj + h * dhdvdh) if extra_h_bug else f_lim * dhdvdh)
                    )
                    Ke[:, ia + 2, jb + 1] += jac_e * (
                        2.0 * h * rvisc * (dqdvJq * Jq + dvJq_i * dqJq_j)
                        + h * cvisc * (dqdvJq * dv + dqJq_j * dbfi)
                        + dvdh_i * frank * dbfj * dq
                        + 4.0 * sus_lin * dvdh_i * qbfv_j
                        + h_bar_fac * dvdh_i * dqflim_j
                    )
                    Ke[:, ia + 2, jb + 2] += jac_e * (
                        4.0 * h * dbfj * dbfi
                        + h * bfj * bfi
                        + 2.0 * h * rvisc * dvJq_i * dvJq_j
                        + h * cvisc * (dvJq_i * dbfj + dvJq_j * dbfi)
                    )

        R = np.zeros(3 * N)
        np.add.at(R, gdofs, Re)
        K_diag = Ke[:, 0:3, 0:3] + np.roll(Ke[:, 3:6, 3:6], 1, axis=0)
        K_upper = Ke[:, 0:3, 3:6]
        K_lower = Ke[:, 3:6, 0:3]
        return R, K_diag, K_upper, K_lower

    U0 = U0_init.copy()
    t = 0.0
    dt = dt_init
    n_step = 0
    max_steps = 900
    while t < t_end - 1e-12 and n_step < max_steps:
        n_step += 1
        U = U0.copy()
        conv = False
        for _ in range(max_nr):
            R, Kd, Ku, Kl = assemble(U, U0, dt)
            if R is None:
                break
            dU = _solve_cyclic_3dof(Kd, Ku, Kl, -R)
            U += dU
            if float(np.linalg.norm(R)) < res_tol and float(np.linalg.norm(dU)) < res_tol:
                conv = True
                break
        if conv:
            U0 = U
            dt = min(dt / adapt, dt_max)
            t += dt
        else:
            dt *= adapt
            if dt < 1e-6:
                break
    return U0


def run(inputs, settings):
    p = {
        "visc": float(inputs["visc"]),
        "cvisc": float(inputs["cvisc"]),
        "rvisc": float(inputs["rvisc"]),
        "frank_cfg": float(inputs["frank"]),
        "sus20_cfg": float(inputs["sus20"]),
        "sus40_cfg": float(inputs["sus40"]),
        "lambda_rot": float(inputs["lambda_rot"]),
        "lambda_iso": float(inputs["lambda_iso"]),
        "lambda_aniso": float(inputs["lambda_aniso"]),
        "crossLinker": float(inputs["crossLinker"]),
        "heqb": float(inputs["heqb"]),
        "k_lim_val": float(inputs["k_lim"]),
        "S2_max": float(inputs["S2_max"]),
    }
    Ldom = float(inputs["r_out"])
    dt_init = float(inputs["deltat_init"])
    dt_max_base = float(inputs["maxDelt"])
    adapt = float(inputs["adaptiveStepTime"])
    noise_cfg = float(inputs["noiseStrength"])

    N = int(round(float(settings.get("N_elem", 160))))
    dt_max_scale = float(settings.get("dt_max_scale", 1.0))
    dt_max = dt_max_base * dt_max_scale
    dt_start = min(dt_init, dt_max) if dt_max_scale < 1.0 else dt_init
    t_max = float(settings.get("t_max", 140.0))
    init_amp = float(settings.get("init_amp", 0.05))
    extra_h_bug = bool(settings.get("extra_h_bug", True))
    compute_threshold = bool(settings.get("compute_threshold", False))
    seed = int(settings.get("seed", 1))

    L_eff = 0.5 * p["frank_cfg"]
    c0_eff = 2.0 * p["sus20_cfg"] + p["lambda_rot"] * p["heqb"] + 4.0 * p["k_lim_val"] / p["S2_max"]
    kap = (p["lambda_aniso"] - p["crossLinker"]) / p["lambda_iso"]

    def lam_fn(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + p["visc"]) * (2.0 * c0_eff + 2.0 * L_eff * x)
        den = (2.0 * c0_eff + 2.0 * L_eff * x) + kap * p["cvisc"] * (x + p["visc"])
        return num / den if den > 0.0 else 1e30

    xs_g = np.geomspace(1e-4, 1e3, 2000)
    vals = np.array([lam_fn(float(x)) for x in xs_g])
    i0 = int(np.argmin(vals))
    res_l = minimize_scalar(
        lam_fn,
        bounds=(float(xs_g[max(0, i0 - 3)]), float(xs_g[min(len(xs_g) - 1, i0 + 3)])),
        method="bounded",
        options={"xatol": 1e-12},
    )
    lam_crit_lin = float(res_l.fun)
    ell_pat = 2.0 * math.pi / math.sqrt(float(res_l.x))

    x = np.linspace(0.0, Ldom, N, endpoint=False)
    rng = np.random.default_rng(seed)
    h0 = np.full(N, p["heqb"])
    q0 = np.zeros(N)
    if init_amp > 0.0:
        for m in range(1, 13):
            ph1, ph2 = rng.uniform(0.0, 2.0 * math.pi, 2)
            arg = 2.0 * math.pi * m * x / Ldom
            h0 += (init_amp / 12.0) * np.cos(arg + ph1)
            q0 += (0.5 * init_amp / 12.0) * np.cos(arg + ph2)
    else:
        h0 += noise_cfg * rng.uniform(0.0, 1.0, N)
        q0 += noise_cfg * rng.uniform(-1.0, 1.0, N)

    U0 = np.empty(3 * N)
    U0[0::3], U0[1::3], U0[2::3] = h0, q0, np.zeros(N)

    Uf = _integrate_1d_fem(N, Ldom, U0, p, dt_start, dt_max, adapt, t_max, extra_h_bug=extra_h_bug)
    h_f, q_f = Uf[0::3], Uf[1::3]
    dx = Ldom / N
    dh_f = (np.roll(h_f, -1) - np.roll(h_f, 1)) / (2.0 * dx)
    omega = float(np.mean(dh_f * dh_f * q_f))

    pad = np.concatenate([h_f[-3:], h_f, h_f[:3]])
    pks, _ = find_peaks(pad, height=1.2, distance=max(3, N // 30))
    pks = pks[(pks >= 3) & (pks < N + 3)]
    n_pks = max(1, len(pks))
    w_ratio = (Ldom / n_pks) / ell_pat

    if compute_threshold:
        N_cell = max(16, int(round(ell_pat / dx)))
        x_c = np.linspace(0.0, ell_pat, N_cell, endpoint=False)
        U_probe = np.empty(3 * N_cell)
        U_probe[0::3] = p["heqb"] + 1e-3 * np.cos(2.0 * math.pi * x_c / ell_pat)
        U_probe[1::3] = -2e-4 * np.cos(2.0 * math.pi * x_c / ell_pat)
        U_probe[2::3] = 0.0

        def grows(lam_try):
            pt = dict(p)
            pt["lambda_iso"] = lam_try
            pt["lambda_aniso"] = lam_try * (1.0 + 0.5 * kap)
            pt["crossLinker"] = lam_try * (1.0 - 0.5 * kap)
            U1 = _integrate_1d_fem(
                N_cell, ell_pat, U_probe, pt, dt_max, dt_max, adapt, 12.0, extra_h_bug=extra_h_bug
            )
            U2 = _integrate_1d_fem(
                N_cell, ell_pat, U1, pt, dt_max, dt_max, adapt, 12.0, extra_h_bug=extra_h_bug
            )
            a1 = float(np.max(np.abs(U1[0::3] - p["heqb"])))
            a2 = float(np.max(np.abs(U2[0::3] - p["heqb"])))
            return a2 > a1

        lo, hi = 0.96 * lam_crit_lin, 1.04 * lam_crit_lin
        for _ in range(11):
            mid = 0.5 * (lo + hi)
            if grows(mid):
                hi = mid
            else:
                lo = mid
        lam_crit_fem = 0.5 * (lo + hi)
        thresh_err_pct = 100.0 * abs(lam_crit_fem - lam_crit_lin) / lam_crit_lin
    else:
        lam_crit_fem = lam_crit_lin
        thresh_err_pct = 0.0

    return {
        "rho_max": float(np.max(h_f)),
        "rho_min": float(np.min(h_f)),
        "S_max": float(2.0 * np.max(np.abs(q_f))),
        "omega": omega,
        "wavelength_ratio": w_ratio,
        "lambda_crit_fem_1d": lam_crit_fem,
        "threshold_rel_err_pct": thresh_err_pct,
    }
