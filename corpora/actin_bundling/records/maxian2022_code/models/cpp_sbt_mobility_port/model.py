"""Faithful Python port of SingleFiberMobilityEvaluator.cpp (lines 122-137, 260-399),
VectorMethods.cpp (lines 228-244), and SingleFiberSaddlePointSolver.cpp (lines 137-161)
in SBT mode (RPYSpecialQuad=RPYDirectQuad=RPYOversample=False) on the Chebyshev kind-2
grid of FibCollocationDiscretization.py (N_x = N_tau + 1)."""

import numpy as np


def _cheb_pts_wts(N, dom, kind):
    a_d, b_d = float(dom[0]), float(dom[1])
    if kind == 1:
        th = np.flipud(2.0 * np.arange(N) + 1.0) * np.pi / (2.0 * N)
        m = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, N, 2) ** 2))
        if N % 2 == 1:
            c = np.concatenate((m, -m[(N + 1) // 2 - 1 : 0 : -1]))
        else:
            c = np.concatenate((m, [0.0], -m[N // 2 - 1 : 0 : -1]))
        v = np.exp(1j * np.arange(N) * np.pi / N)
        w = np.real(np.fft.ifft(c * v))
    elif kind == 2:
        th = np.flipud(np.pi * np.arange(N)) / (N - 1)
        c = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, N, 2) ** 2))
        c = np.concatenate((c, c[N // 2 - 1 : 0 : -1]))
        w = np.real(np.fft.ifft(c))
        w[0] /= 2.0
        w = np.concatenate((w, [w[0]]))
    else:
        raise ValueError(f"Unsupported kind: {kind}")
    x = 0.5 * (b_d - a_d) * (np.cos(th) + 1.0) + a_d
    w = 0.5 * (b_d - a_d) * w
    return x, w, th


def _resamp_mat(Ntarg, Nsrc, kind_t, kind_s):
    _, _, th_s = _cheb_pts_wts(Nsrc, [-1.0, 1.0], kind_s)
    _, _, th_t = _cheb_pts_wts(Ntarg, [-1.0, 1.0], kind_t)
    C_s = np.cos(np.outer(th_s, np.arange(Nsrc)))
    C_t = np.cos(np.outer(th_t, np.arange(Nsrc)))
    return C_t @ np.linalg.inv(C_s)


def _diff_mat(N, dom, kind):
    _, _, th = _cheb_pts_wts(N, dom, kind)
    CtoV = np.cos(np.outer(th, np.arange(N)))
    VtoC = np.linalg.inv(CtoV)
    Dcos = np.zeros_like(VtoC)
    Dcos[N - 2, :] = 2.0 * (N - 1) * VtoC[N - 1, :]
    for k in range(N - 3, -1, -1):
        Dcos[k, :] = Dcos[k + 2, :] + 2.0 * (k + 1) * VtoC[k + 1, :]
    Dcos[0, :] *= 0.5
    return (CtoV @ Dcos) / (0.5 * (dom[1] - dom[0]))


def _regularize_s(sX, delta, L):
    # FibCollocationDiscretization.py:168-176
    sNew = 0.5 * L * np.ones_like(sX)
    if delta < 0.5:
        x = 2.0 * sX / L - 1.0
        regwt = np.tanh((x + 1.0) / delta) - np.tanh((x - 1.0) / delta) - 1.0
        sNew = sX.copy()
        mask = sX < 0.5 * L
        sNew[mask] = regwt[mask] * sX[mask] + (1.0 - regwt[mask] ** 2) * delta * L / 2.0
        sNew[sX > 0.5 * L] = L - np.flip(sNew[mask])
    return sNew


def run(inputs, settings):
    eps = float(inputs["eps"])
    L = float(inputs["L"])
    # Unit note: 1 Pa*s = 1 pN*s/um^2 identically!
    mu = float(inputs["mu"])
    delta = float(inputs["delta"])
    eig_thres = float(inputs["eig_val_thres"])
    N_tau = int(round(float(settings.get("N_tau", 16))))
    Nx = N_tau + 1

    # FibCollocationDiscretization.py:10, 105
    a_rpy = np.exp(1.5) / 4.0 * eps * L

    # Energy discretization uses kind-2 grid for X (FibCollocationDiscretization.py:274-277)
    sX, _, thX = _cheb_pts_wts(Nx, [0.0, L], 2)
    DX = _diff_mat(Nx, [0.0, L], 2)
    sbar = _regularize_s(sX, delta, L)

    # SingleFiberMobilityEvaluator.cpp:260-299 (initLocalDragCoeffsRPY & updateLocalDragCoeffsRPY, nXs=1)
    loc_s = np.where(
        sbar < 2.0 * a_rpy,
        np.log((L - sbar) / (2.0 * a_rpy)),
        np.where(
            sbar > L - 2.0 * a_rpy,
            np.log(sbar / (2.0 * a_rpy)),
            np.log(sbar * (L - sbar) / (4.0 * a_rpy * a_rpy)),
        ),
    )
    loc_d = np.where(
        sbar < 2.0 * a_rpy,
        1.0 / (8.0 * a_rpy * a_rpy) - 1.0 / (2.0 * (L - sbar) ** 2),
        np.where(
            sbar > L - 2.0 * a_rpy,
            1.0 / (8.0 * a_rpy * a_rpy) - 1.0 / (2.0 * sbar ** 2),
            1.0 / (4.0 * a_rpy * a_rpy) - 1.0 / (2.0 * sbar ** 2) - 1.0 / (2.0 * (L - sbar) ** 2),
        ),
    )
    loc_I = np.where(
        sbar < 2.0 * a_rpy,
        (128.0 * a_rpy ** 2 - 36.0 * a_rpy ** 2 + 64.0 * a_rpy * sbar - 9.0 * sbar ** 2) / (48.0 * a_rpy ** 2),
        np.where(
            sbar > L - 2.0 * a_rpy,
            (128.0 * a_rpy ** 2 - 36.0 * a_rpy ** 2 + 64.0 * a_rpy * (L - sbar) - 9.0 * (L - sbar) ** 2) / (48.0 * a_rpy ** 2),
            (32.0 - 9.0) / 6.0,
        ),
    )
    loc_T = np.where(
        sbar < 2.0 * a_rpy,
        0.25 + sbar ** 2 / (16.0 * a_rpy ** 2),
        np.where(
            sbar > L - 2.0 * a_rpy,
            0.25 + (L - sbar) ** 2 / (16.0 * a_rpy ** 2),
            0.5,
        ),
    )

    visc_inv = 1.0 / (8.0 * np.pi * mu)
    # SingleFiberMobilityEvaluator.cpp:324-326 (tau = [1,0,0])
    Mloc_par = visc_inv * (2.0 * loc_s - (4.0 * a_rpy ** 2 / 3.0) * loc_d + loc_I + loc_T)
    Mloc_perp = visc_inv * (loc_s + (2.0 * a_rpy ** 2 / 3.0) * loc_d + loc_I)

    # FibCollocationDiscretization.py:599-635 & SingleFiberMobilityEvaluator.cpp:332-399
    sscale = -1.0 + 2.0 * sX / L
    CtoV = np.cos(np.outer(thX, np.arange(Nx)))
    AllQs = np.zeros((Nx, Nx))
    for iPt in range(Nx):
        eta = sscale[iPt]
        q = np.zeros(Nx)
        if eta > -1.0:
            n_lo, w_lo, _ = _cheb_pts_wts(200, [-1.0, eta], 1)
            ac = np.arccos(np.clip(n_lo, -1.0, 1.0))
            for kk in range(Nx):
                q[kk] += np.dot(w_lo, -np.cos(kk * ac))
        if eta < 1.0:
            n_hi, w_hi, _ = _cheb_pts_wts(200, [eta, 1.0], 1)
            ac = np.arccos(np.clip(n_hi, -1.0, 1.0))
            for kk in range(Nx):
                q[kk] += np.dot(w_hi, np.cos(kk * ac))
        AllQs[iPt, :] = q
    FPMatrix = visc_inv * 0.5 * L * np.linalg.solve(CtoV.T, AllQs.T).T

    FP_scalar = np.zeros((Nx, Nx))
    for iPt in range(Nx):
        for jPt in range(Nx):
            if iPt == jPt:
                FP_scalar[iPt, :] += FPMatrix[iPt, iPt] * DX[iPt, :]
            else:
                ds = sX[jPt] - sX[iPt]
                FP_scalar[iPt, jPt] += (1.0 / ds) * FPMatrix[iPt, jPt]
                FP_scalar[iPt, iPt] -= (1.0 / ds) * FPMatrix[iPt, jPt]

    M_par = np.diag(Mloc_par) + 2.0 * FP_scalar
    M_perp = np.diag(Mloc_perp) + 1.0 * FP_scalar

    # FibCollocationDiscretization.py:520-524 -> W_tilde = R^T W_{2Nx} R
    _, w2, _ = _cheb_pts_wts(2 * Nx, [0.0, L], 2)
    R2 = _resamp_mat(2 * Nx, Nx, 2, 2)
    Wt = R2.T @ np.diag(w2) @ R2
    WtInv = np.linalg.inv(Wt)

    # SingleFiberMobilityEvaluator.cpp:134-135 & VectorMethods.cpp:228-244
    def _sym_and_clamp(M):
        MF = M @ WtInv
        MF_sym = 0.5 * (MF + MF.T)
        ev, evec = np.linalg.eigh(MF_sym)
        raw_min = float(ev[0])
        ev_clamped = np.maximum(ev, eig_thres)
        return evec @ np.diag(ev_clamped) @ evec.T, raw_min

    MF_par, min_par = _sym_and_clamp(M_par)
    MF_perp, min_perp = _sym_and_clamp(M_perp)

    ones = np.ones(Nx)
    sc = sX - 0.5 * L
    N_xx = 1.0 / float(ones @ np.linalg.solve(MF_par, ones))
    N_yy = 1.0 / float(ones @ np.linalg.solve(MF_perp, ones))
    N_rot = 1.0 / float(sc @ np.linalg.solve(MF_perp, sc))

    return {
        "alpha": float(mu * L * N_yy),
        "beta": float(mu * L * (N_xx - N_yy)),
        "gamma": float(mu * (L ** 3) * N_rot),
        "min_raw_eig": float(min(min_par, min_perp)),
    }
