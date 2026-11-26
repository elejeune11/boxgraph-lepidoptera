"""Faithful Python port of the SM Sec. C local-drag rigid-rod mobility (RegularizeS.m +
AsymptoticLocalMobility.m, local drag only) coupled to CrossLinkForceEvaluator.cpp::calcCLForces
(smooth_force=True) or calcCLForcesEnergy (smooth_force=False) to compute the exact linearized
explicit-Euler stability limit dt_max for two parallel rigid rods linked at s = L."""

import numpy as np


def _cheb_kind2(Nx, L):
    th = np.flipud(np.pi * np.arange(Nx)) / (Nx - 1)
    c = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, Nx, 2) ** 2))
    c = np.concatenate((c, c[Nx // 2 - 1 : 0 : -1]))
    w = np.real(np.fft.ifft(c))
    w[0] /= 2.0
    w = np.concatenate((w, [w[0]]))
    s = 0.5 * L * (np.cos(th) + 1.0)
    w = 0.5 * L * w
    return s, w


def _regularize_s(s0, delta, L):
    sNew = 0.5 * L * np.ones_like(s0)
    if 0.0 < delta < 0.5:
        x = 2.0 * s0 / L - 1.0
        regwt = np.tanh((x + 1.0) / delta) - np.tanh((x - 1.0) / delta) - 1.0
        sNew = s0.copy()
        mask = s0 < 0.5 * L
        sNew[mask] = regwt[mask] * s0[mask] + (1.0 - regwt[mask] ** 2) * delta * L / 2.0
        sNew[s0 > 0.5 * L] = L - np.flip(sNew[mask])
    return sNew


def run(inputs, settings):
    L = float(inputs["L"])
    eps = float(inputs["eps"])
    mu = float(inputs["mu"])
    K_c = float(inputs["K_c"])
    delta = float(inputs["delta"])
    sigma_rel = float(inputs["sigma_rel"])
    dt_small = float(inputs["dt_small"])
    dt_large = float(inputs["dt_large"])

    N = int(round(float(settings.get("N", 50))))
    smooth_force = bool(settings.get("smooth_force", True))
    Nx = N + 1

    s, w = _cheb_kind2(Nx, L)
    sbar = _regularize_s(s, delta, L)
    a = eps * L
    c_val = np.log(4.0 * sbar * (L - sbar) / (a * a))
    Mloc_perp = (c_val + 1.0) / (8.0 * np.pi * mu)

    sc = s - 0.5 * L
    N_perp = 1.0 / float(np.dot(w, 1.0 / Mloc_perp))
    N_rot = 1.0 / float(np.dot(w, (sc ** 2) / Mloc_perp))

    if smooth_force:
        # CrossLinkForceEvaluator.cpp:96-140, 273-280 at s1* = s2* = L
        sig = sigma_rel * L
        dh = np.exp(-((s - L) ** 2) / (2.0 * sig * sig)) / np.sqrt(2.0 * np.pi * sig * sig)
        dh_norm = dh / float(np.dot(w, dh))
        arm = float(np.dot(w, sc * dh_norm))
        # Exact Jacobian on antisymmetric mode (y1 = -y2, th1 = -th2) for rigid rods:
        # F_1 = -2 K_c (y + arm*th), T_1 = arm * F_1 (since t0 x t0 = 0 kills m2 term)
        m_end = N_perp + (arm ** 2) * N_rot
    else:
        arm = 0.5 * L
        m_end = N_perp + (arm ** 2) * N_rot

    # Stability of explicit Euler: 1 - 2 * K_c * m_end * dt >= -1 => dt_max = 1 / (K_c * m_end)
    dt_max = 1.0 / (K_c * m_end)

    return {
        "alpha_ld": float(mu * L * N_perp),
        "gamma_ld": float(mu * (L ** 3) * N_rot),
        "effective_arm": float(arm),
        "dt_max": float(dt_max),
        "dt_small_pct": float(100.0 * dt_small / dt_max),
        "dt_large_pct": float(100.0 * dt_large / dt_max),
    }
